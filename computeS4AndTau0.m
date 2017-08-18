function [ S4, tau0 ] = computeS4AndTau0( zkhist, tkhist )
%COMPUTES4ANDTAU0 Compute scintillation index and decorrelation time
%   [S4, TAU0] = COMPUTES4ANDTAU0(Z, T) computes the scintillation index
%   S4 and the decorrelation time TAU0 corresponding to the input complex
%   channel response function time history Z
%
%   INPUTS
%
%   zkhist      Nt-by-1 vector containing the normalized complex scintillation
%               time history in the form of averages over Ts with sampling
%               interval Ts.  zkhist(kp1) is the average over tk to tkp1.
%
%   tkhist      Nt-by-1 vector of time points corresponding to zkhist.
%
%
%   OUTPUTS
%
%   S4          Intensity scintillation index of the scintillation time history
%               in zkhist, equal to the mean-normalized standard deviation of
%               the intensity abs(zkhist).^2.
%
%   tau0        The decorrelation time of the scintillation time history in
%               zkhist, in seconds.
%
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+

I = abs(zkhist).^2;
S4 = std(I) / mean(I);

% Find time-varying component of Z
xsi = zkhist - mean(zkhist);

% Autocorrelate
[R, ii] = ccorr(xsi, xsi);

% Only consider positive time values
R = R(ii >= 0);

% Find time where peak value has decayed by 1/e
R0 = R(1);
[~, idx] = min(abs(R0*exp(-1) - R));
tau0 = tkhist(idx);

end

