function [ Ad, Bd, Cd, Dd, BnAct ] = configureloopfilter( Bn, Ta, looporder )
%CONFIGURELOOPFILTER Configure a discrete-time loop filter for a feedback
%tracking loop.
%   [AD,BD,CD,DD,BNACT] = CONFIGURELOOPFILTER(BN,TA,LOOPORDER) configures
%   a discrete-time loop filter of order LOOPORDER with sampling interval
%   TA and target loop noise bandwidth BN. AD, BD, CD, and DD make up the
%   state-space representation of the discrete system, and BNACT is the actual 
%   loop noise bandwidth.

switch looporder
    case 1
        K = 4*Bn;
        Ds = K * tf(1, 1);
    case 2
        K = 8*Bn/3;
        a = K/2;
        Ds = K * tf([1 a], [1 0]);
    case 3
        a = 1.2 * Bn;
        b = a^2 / 2;
        K = 2*a;
        Ds = K * tf([1 a b], [1 0 0]);
    otherwise
        error('Loop order must be 1, 2, or 3');
end

% Convert the loop filter to a discrete-time state-space model
Dz = c2d(Ds, Ta, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(Dz);

NCO = tf(Ta, [1 -1], Ta);
PD = tf([1 1], [2 0], Ta);
Fz = PD * Dz * NCO;
Hz = Fz/(1 + Fz);

% Calculate BnAct
wAlias = pi/Ta;
wVec = (0:10000)' * (wAlias/10000);
[magVec, ~] = bode(Hz, wVec);
magVec = magVec(:);
BnAct = sum(magVec.^2)*mean(diff(wVec))/(2*pi);

end

