function [ Sp, Se, Sl ] = gnsscorrelate( x, t0, fIF, Fs, Ta, ts, fD, thetaHat, teml, code )
%GNSSCORRELATE Complex correlation with a bandpass GNSS signal.
%   [SP,SE,SL] = GNSSCORRELATE(X,T0,FIF,FS,TA,TS,FD,THETAHAT,TEML,PRN)
%   calculates the prompt, early, and late complex accumulations over one
%   accumulation interval. X is one accumulation interval of the bandpass data 
%   to be correlated against, T0 is the receiver time corresponding to the
%   first sample in X, FIF is the intermediate (center) bandpass frequency,
%   FS is the sample rate of X, TA is the accumulation interval, TS is the start
%   time of the corresponding C/A code (in receiver time), FD is the current 
%   Doppler estimate, THETAHAT is the initial phase offset estimate for the 
%   accumulation, TEML is the early-minus-late time (in seconds), and PRN is the
%   unique 1023 length C/A code of the satellite beinig tracked.
%
%   SP, SE, and SL are the prompt, early, and late accumulations,
%   respectively.

Nc = 1023;
Pc = Ta / round(Ta / 1e-3);
Nk = length(x);
Tc =  Pc / Nc;

numCodes = round((Nk / Fs) / Pc);
code = repmat(code, [numCodes 1]);

codePrompt = oversampleSpreadingCode(...
    code, ...
    (1 / Fs) / Tc, ...
    (t0 - ts) / Tc, ...
    Nk, ...
    numCodes * Nc);

codeEarly = oversampleSpreadingCode(...
    code, ...
    (1 / Fs) / Tc, ...
    (t0 - (ts - teml)) / Tc, ...
    Nk, ...
    numCodes * Nc);

codeLate = oversampleSpreadingCode(...
    code, ...
    (1 / Fs) / Tc, ...
    (t0 - (ts + teml)) / Tc, ...
    Nk, ...
    numCodes * Nc);

t = t0 + (0:Nk-1)' * 1/Fs;
loc = exp(-1j * (2*pi*(fIF*t + fD*(t - t0)) + thetaHat));
Sp = sum(x .* loc .* codePrompt);
Se = sum(x .* loc .* codeEarly);
Sl = sum(x .* loc .* codeLate);

end

