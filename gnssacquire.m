=function [ fD, ts, CN0, sigma, Z ] = gnssacquire( x, fc, fs, txId, Ta, fspan, N, showPlot )
%GNSSACQUIRE Acquire GNSS signals from bandpass data.
%   [FD,TS] = GNSSACQUIRE(X,FC,FS,TXID,TA) searches for
%   the optimal values for FD (Hz) and TS (sec) that maximize the correlation 
%   between the bandpass signal X with center frequency FC and the PRN sequence 
%   for satellite TXID. If TXID is an array, the acquisition is repeated for 
%   each element of TXID. FS is the sampling rate of the data in X and TA is the
%   acquisition interval over which to perform the search. If X is baseband
%   data, set FC=0.
%
%   [FD,TS] = GNSSACQUIRE(X,FC,FX,TXID,TA,FSPAN) specifies
%
%   [FD,TS,CN0,S] = GNSSACQUIRE(...) also computes the carrier-to-noise
%   ratio CNO (in dB-Hz) and the standard deviation of the noise in the
%   signal S.
%
%   [FD,TS,...] = GNSSACQUIRE(...,1) or GNSSACQUIRE(...) (with no output 
%   arguments) plots the squared magnitude of the correlation along the FD and 
%   TS dimensions.
%
%   The signal X is assumed to be low-side mixed. If this is not the case, the
%   Doppler estimate FD must be multiplied by -1.

% Set defaults if arguments missing
if nargin < 6
    fspan = 7000;
    N = 1;
    showPlot = 0;
elseif nargin < 7
    N = 1;
    showPlot = 0;
elseif nargin < 8
    showPlot = 0;
end

% Define constants
CODE_PERIOD = 1e-3;
CODE_CHIPS = 1023;
CODE_CHIP_INTERVAL = CODE_PERIOD / CODE_CHIPS;

nCode = floor(fs * CODE_PERIOD);               % Samples per code
nAcq = floor(fs * Ta);                         % Samples per acquisition
nCodesPerAcq = floor(nAcq / nCode);            % Codes per acquisition period
nFft = 2^nextpow2(nAcq);

xAcq = x(1:N*nFft);
tVec = (1/fs) * (0:nFft-1)';

% Frequency search vector
if length(fspan) == 1
    dfD = 1/(4*Ta);
    fDSearchVec = -fspan:dfD:fspan;
elseif length(fspan) == 2
    dfD = 1/(4*Ta);
    fDSearchVec = fspan(1):dfD:fspan(2);        
else
    dfD = mean(diff(fspan));
    fDSearchVec = fspan;
end

% Calculate all values at once instead of in every loop for efficiency
if fc > 0
    localMat = exp(-1j * (2*pi*tVec*(fc + fDSearchVec)));
elseif fc == 0
    localMat = exp(-1j * (tVec * 2*pi*fDSearchVec));
else
    error('Center frequency must be greater than or equal to 0.');
end

fD = zeros(size(txId));
ts = zeros(size(txId));
CN0 = zeros(size(txId));
sigma = zeros(size(txId));
Z = zeros(size(txId));

for m = 1:length(txId)
    prn = oversampleSpreadingCode(...
            repmat(generatePrnSeq(txId(m)), [nCodesPerAcq 1]), ...
            (1 / fs) / CODE_CHIP_INTERVAL, ...
            0, ...
            nAcq, ...
            nCodesPerAcq * CODE_CHIPS);

    prnFft = conj(fft(prn, nFft));

    for n = 1:N
        M = zeros(nCode, length(fDSearchVec));
        for i = 1:length(fDSearchVec)
            iidx = (n-1)*nAcq + (1:nFft);
            xFft = fft(xAcq(iidx) .* localMat(:, i));
            yFft = xFft .* prnFft;
            y = ifft(yFft);
            M(:, i) = abs(y(1:nCode)).^2;
        end

        S = max(max(M));
        
        if n == 1
            [maxRow, maxCol] = find(M == S);
            fDVec = M(maxRow, :);
            tsVec = M(:, maxCol);
            
            if nargout > 0
                fD(m) = fDSearchVec(maxCol);
                jk = mod(maxRow - 1, nCode);
                ts(m) = jk/fs;

                if nargout > 2
                    % To calculate C/N0 we must first find sigmaIQ. To do that, we must cut out
                    % values close to the correct Doppler frequency and time offset
                    idxf = 2/(Ta * dfD);
                    idxt = ceil(CODE_CHIP_INTERVAL * fs);
                    iif = max(maxCol-idxf, 1):min(maxCol+idxf, length(fDSearchVec));
                    iit = max(maxRow-idxt, 1):min(maxRow+idxt, nCode);
                    M(iit, maxCol) = 0;
                    M(maxRow, iif) = 0;

                    twosigmaIQ2 = sum(sum(M))/sum(sum(M ~= 0));
                    CN0(m) = 10*log10((S - twosigmaIQ2)/(twosigmaIQ2*Ta));
                    sigma(m) = sqrt(twosigmaIQ2/2);
                end
            end
        end
        
        Z(m) = Z(m) + S;
    end
    
    if showPlot || nargout == 0
        figure(txId(m));
        clf;
        subplot(1, 2, 1), ...
            hold on, ...
            plot(fDSearchVec, fDVec), ...
%             plot([fDSearchVec(1) fDSearchVec(end)], [1 1] * 1e7, 'r'), ...
            xlabel('f_D (Hz)'), axis tight;
        subplot(1, 2, 2), ...
            hold on, ...
            plot(tsVec), ...
%             plot([0 length(tsVec)-1], [1 1] * 1e7, 'r'), ...
            xlabel('Code offset (samples)'), axis tight;
    end
end
