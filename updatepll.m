function [ xkp1, vk ] = updatepll( s )
%UPDATEPLL Perform a single update step of a phase tracking loop with an
%arctangent phase detector
%
%   INPUTS
%
%   s ------------- A structure with the following fields:
%
%       Ip --------- The in-phase prompt accumulation over the interval from
%                    tkm1 to tk.
%
%       Qp --------- The quadrature prompt accumulation over the interval from
%                    tkm1 to tk.
%
%       xk --------- The phase tracking loop filter's state at time tk. The
%                    dimension of xk is N-1, where N is the order of the loop's
%                    closed-loop transfer function.
%
%       Ad,Bd,Cd,Dd -- The loop filter's state-space model.
%
%   OUTPUTS
%
%   xkp1 -------- The loop filter's state at time tkp1. The dimension of xkp1
%                 is N-1, where N is the order of the loop’s closed-loop
%                 transfer function.
%
%   vk ---------- The Doppler frequency shift that will be used to drive the
%                 receiver's carrier-tracking numerically controlled
%                 oscillator during the time interval from tk to tkp1, in
%                 rad/sec.
%
%+-----------------------------------------------------------------------------+
%   References:
%
%
%+=============================================================================+

ek = atan(s.Qp / s.Ip);

xkp1 = s.Ad * s.xk + s.Bd * ek;
vk = s.Cd * s.xk + s.Dd * ek;

end

