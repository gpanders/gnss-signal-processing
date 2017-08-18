function [ x ] = iq2if( I, Q, Tl, fIF )
%IQ2IF Convert baseband I and Q samples to intermediate frequency samples
%   Let xl(m*Tl) = I(m*Tl) + j*Q(m*Tl) be a discrete-time baseband
%   representation of a bandpass signal. This function converts xl(n) to a
%   discrete-time bandpass signal x(n) = I(n*T)*cos(2*pi*fIF*n*T) -
%   Q(n*T)*sin(2*pi*fIF*n*T) centered at the user-specified intermediate
%   frequency fIF, where T = Tl/2.
%
%   INPUTS
%
%   I ------------- N-by-1 vector of in-phase baseband samples.
%
%   Q ------------- N-by-1 vector of quadrature baseband samples.
%
%   Tl ------------ Sampling interval of baseband samples (complex sampling
%                   interval), in seconds.
%
%   fIF ----------- Intermediate frequency to which the baseband samples will
%                   be up-converted, in Hz
%
%
%   OUTPUTS
%
%   x ------------- 2*N-by-1 vector of intermediate frequency samples with
%                   sampling interval T = Tl/2
%
%+------------------------------------------------------------------------------+
% References:
%
%
%+==============================================================================+

if length(I) ~= length(Q)
    error('I and Q must be same length');
end

N = length(I);
t = (0:2*N-1)' * Tl/2;
xl = sqrt(2) * exp(1j*2*pi*fIF*t);

S = interp(I + 1i*Q, 2);
x = real(S .* xl);

end

