close all; clearvars -except Y tVec;

Fs = 40e6/7;         % Sampling frequency (Hz)
Tfull = 60;         % Time interval of data to load
N = Fs * Tfull;        
N = floor(N/16) * 16;  % Number of data samples to load

fid = fopen('mysteryData3.bin', 'r','l');
[Y, count] = binloadSamples(fid, N, 'dual');
Y = Y(:,1);
if(count ~= N)
  error('Insufficient data');
end

% dfDataHead70 contains 70 seconds of IF data from dfDataHead.bin
% load('dfDataHead70.mat');

% mysteryData3 contains 60 seconds of IF data from mysteryData3.bin
% load('mysteryData3.mat');

% mysteryData2 contains 40 seconds of IF data from mysteryData2.bin
% load('mysteryData2.mat');

% Define constants
fL1 = 154 * 10.23e6;
fIF = 1610476.19047612;
Ts = 1 / Fs;
Nc = 1023;                          % number of chips per code period
Tc = 1e-3/Nc;                       % code chip interval
Pc = Nc * Tc;                       % code period
Ta = 1e-3;                          % acquisition interval
t0 = 0;                             % initial time
tf = Tfull;
teml = 0.5 * Tc;                    % early-minus-late time
tVec = (0:N-1)' * Ts;               % time vector
sMix = 1;                          % -1 for high-side mixed, 1 for low-side mixed

showPlots = 1;

numAcq = floor((tf - t0)/Ta);  % number of acquisitions to perform

txIdVec = [1 7 8 11 28 30]';

% Calculate indices for start and end times, and update times accordingly
% so they correspond to an actual time in `tVec`
t0idx = max(ceil(t0 * Fs), 1);
t0 = tVec(t0idx);

tfidx = min(floor(tf * Fs), N);
tf = tVec(tfidx);

% Preallocate time history variables
tsVec = zeros(numAcq+1, length(txIdVec));
fDVec = zeros(numAcq+1, length(txIdVec));
thetaHatVec = zeros(numAcq+1, length(txIdVec));
IVec = zeros(numAcq, length(txIdVec));
QVec = zeros(numAcq, length(txIdVec));
CN0Vec = zeros(numAcq, length(txIdVec));

% Initialize DLL
dll.Bn = 0.1;

% Initialize PLL
pll.Bn = 10;
[pll.Ad, pll.Bd, pll.Cd, pll.Dd, ~] = configureloopfilter(pll.Bn, Ta, 3);

% Main loop
for m = 1:length(txIdVec)
    txId = txIdVec(m);
    prn = generatePrnSeq(txId);
    
    % Perform coarse acquisition
    [fDcoarse, ~, CN0, sigmaIQ, Z] = gnssacquire(Y(t0idx:tfidx), fIF, Fs, txId, Ta, [-40e3 -10e3]);

    % Perform fine Acquisition
    [fD, ts] = gnssacquire(Y(t0idx:tfidx), fIF, Fs, txId, 10e-3, fDcoarse + (-250:2:250));
    
    % If straddling a data bit, perform a 2nd 10ms acquisition
    % Not sure how to check for this other than manually examining the
    % plots generated by GNSSACQUIRE...
    if ismember(txId, [7 30])
        idxt02 = ceil((t0 + 10e-3) * Fs);
        fD = gnssacquire(Y(idxt02:tfidx), fIF, Fs, txId, 10e-3, fDcoarse + (-250:2:250));
    end

    % Adjust Doppler estimate for mixing type
    fD = sMix * fD;

    % Use acquisition values as initial values
    tsVec(1, m) = t0 + ts;
    fDVec(1, m) = fD;
    
    % Prime the loop filter
    pll.xk = [pll.Cd; pll.Cd*pll.Ad] \ ones(2, 1) * 2*pi*(sMix*fD);
    
    % Configure DLL
    dll.sigmaIQ = sigmaIQ;
    dll.IsqQsqAvg = Z;
    
    % Adjust acquisition interval for Doppler effect
    Tak = Ta/(1 + (sMix*fD)/fL1);
    
    % Find first sample after first code start time
    jk = ceil(tsVec(1, m) * Fs);
    % First sample after second code start time
    jkp1 = jk + ceil(Fs*Tak);
    
    % First estimates of theta and ts
    thetaHatVec(2, m) = thetaHatVec(1, m) + (tVec(jkp1) - tVec(jk))*(sMix * 2*pi*fD);
    tsVec(2, m) = tsVec(1, m) + (1 - (sMix*fD)/fL1)*Pc;
    
    % Start tracking
    for k = 2:numAcq+1
        thetaHat = thetaHatVec(k-1, m);
        fD = fDVec(k-1, m);
        ts = tsVec(k-1, m);
        
        % Perform accumulation
        [Sp, Se, Sl] = gnsscorrelate(Y(jk:jkp1-1), tVec(jk), fIF, Fs, Tak, ts, sMix*fD, thetaHat, teml, prn);

        IVec(k-1, m) = real(Sp);
        QVec(k-1, m) = imag(Sp);

        % Feed values to PLL and update
        pll.Ip = IVec(k-1, m);
        pll.Qp = QVec(k-1, m);

        [xkp1, vk] = updatepll(pll);
        
        pll.xk = xkp1;

        % Update DLL
        dll.Tc = (1 - (sMix*fD)/fL1) * Tc;
        dll.vp = sMix * vk / (2*pi*fL1);
        dll.IsqQsqAvg = ((k-1)*dll.IsqQsqAvg + (abs(Sp)^2))/k;
        dll.Ip = real(Sp);
        dll.Qp = imag(Sp);
        dll.Ie = real(Se);
        dll.Qe = imag(Se);
        dll.Il = real(Sl);
        dll.Ql = imag(Sl);

        vTotal = updatedll(dll);

        CN0Vec(k-1, m) = 10*log10((dll.IsqQsqAvg - 2*sigmaIQ^2)/(2*sigmaIQ^2 * Ta));
        fDVec(k, m) = sMix * vk/(2*pi);

        if k <= numAcq
            thetaHatVec(k+1, m) = thetaHatVec(k, m) + (tVec(jkp1) - tVec(jk))*vk;
            tsVec(k+1, m) = tsVec(k, m) + (1 - vTotal)*Pc;
        end      

        % Calculate accumulation interval for next accumulation
        Tak = Ta/(1 + sMix*fDVec(k, m)/fL1);
        
        % Update indices
        jk = jkp1;
        jkp1 = jk + ceil(Fs*Tak);
        
        % Make sure there is still more data to go
        if jkp1 > tfidx
            % We assumed there would be more acquisitions than there actually
            % was, so just quit
            break;
        end   
    end
    
    if k < numAcq+1
        % We stopped short of how many we thought there would be. Adjust
        % `numAcq' for future iterations (i.e. next txId)
        numAcq = k-1;
        
        % Delete unused cells
        tsVec(k+2:end, :) = [];
        thetaHatVec(k+2:end, :) = [];
        fDVec(k+1:end, :) = [];
        IVec(k:end, :) = [];
        QVec(k:end, :) = [];
        CN0Vec(k:end, :) = [];
    end
    
    if showPlots
        % Plot I's and Q's over time
        figure((m-1)*4+1);
        clf;
        hold on;
        plot(t0 + (0:numAcq-1) * Ta, IVec(:, m), 'g');
        plot(t0 + (0:numAcq-1) * Ta, QVec(:, m), 'b');
        legend('I_p', 'Q_p');
        title(['Satellite ' num2str(txId)]);
        xlim(t0 + [0 numAcq-1]*Ta);
        ylim([-2500 2500]);
        xlabel('Time (s)');
        ylabel('Power');

        % Plot I's and Q's on complex plot with time as color bar
        figure((m-1)*4+2);
        clf;
        scatter(IVec(:, m), QVec(:, m), 10, linspace(0, 255, numAcq), 'filled');
        grid on;
        axis equal;
        xlabel('I');
        ylabel('Q');
        title(['Satellite ' num2str(txId)]);
        h = colorbar;
        tTicks = (0:10:Tfull)';
        h.Ticks = linspace(0, 255, length(tTicks));
        h.TickLabels = cellstr(num2str(tTicks));
        h.Label.String = 'Time (s)';

        % Plot fD estimate over time
        figure((m-1)*4+3);
        clf;
        plot(t0 + (0:numAcq) * Ta, fDVec(:, m));
        xlabel('Time (s)');
        ylabel('f_D (Hz)');
        title(['Satellite ' num2str(txId)]);

        % Plot calculated C/N0 over time
        figure((m-1)*4+4);
        clf;
        plot(t0 + (0:numAcq-1) * Ta, CN0Vec(:, m));
        xlabel('Time (s)');
        ylabel('C/N0 (dB-Hz)');
        ylim([37 51]);
        title(['Satellite ' num2str(txId)]);
    end
end