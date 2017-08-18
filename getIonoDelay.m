function [ delTauG ] = getIonoDelay( ionodata, fc, rRx, rSv, tGPS, model )
% GETIONODELAY Return a model-based estimate of the ionospheric delay
% experienced by a transionospheric GNSS signal as it propagates from a 
% GNSS SV to the antenna of a terrestrial GNSS receiver.
%
% INPUTS
%
% ionodata ------- Structure containing a parameterization of the
%                  ionosphere that is valid at time tGPS.  The structure is
%                  defined differently depending on what ionospheric model
%                  is selected:
%
%                  broadcast --- For the broadcast (Klobuchar) model, ionodata
%                  is a structure containing the following fields:
%
%                     alpha0 ... alpha3 -- power series expansion coefficients
%                     for amplitude of ionospheric TEC
%
%                     beta0 .. beta3 -- power series expansion coefficients
%                     for period of ionospheric plasma density cycle
%
%
%                  Other models TBD ...
%
% fc ------------- Carrier frequency of the GNSS signal, in Hz.
%
% rRx ------------ A 3-by-1 vector representing the receiver antenna position
%                  at the time of receipt of the signal, expressed in meters
%                  in the ECEF reference frame.
%
% rSv ------------ A 3-by-1 vector representing the space vehicle antenna
%                  position at the time of transmission of the signal,
%                  expressed in meters in the ECEF reference frame.
%
% tGPS ----------- A structure containing the true GPS time of receipt of
%                  the signal.  The structure has the following fields:
%
%                  week -- unambiguous GPS week number
%
%                  seconds -- seconds (including fractional seconds) of the
%                  GPS week
%
% model ---------- A string identifying the model to be used in the
%                  computation of the ionospheric delay:
%
%                  broadcast --- The broadcast (Klobuchar) model.
%
%                  Other models TBD ...
%
% OUTPUTS
%
% delTauG -------- Modeled scalar excess group ionospheric delay experienced
%                  by the transionospheric GNSS signal, in seconds.
%
%+------------------------------------------------------------------------------+
% References: For the broadcast (Klobuchar) model, see IS-GPS-200F
% pp. 128-130.
%
%+==============================================================================+

if strcmpi(model, 'broadcast') || strcmpi(model, 'klobuchar')
    alpha = [ionodata.alpha0 ionodata.alpha1 ionodata.alpha2 ionodata.alpha3];
    beta = [ionodata.beta0 ionodata.beta1 ionodata.beta2 ionodata.beta3];
    [E, A] = satelaz(rSv, rRx);
    [lat, lon, ~] = ecef2lla(rRx);

    % Convert to semi-circles
    E = E/pi;
    A = A/pi;
    lat = lat/pi;
    lon = lon/pi;

    psi = 0.0137/(E + 0.11) - 0.022;
    phi_i = lat + psi*cos(A*pi);
    if abs(phi_i) > 0.416
        phi_i = sign(phi_i) * 0.416;
    end

    lambda_i = lon + psi*sin(A*pi)/cos(phi_i*pi);
    phi_m = phi_i + 0.064 * cos((lambda_i - 1.617)*pi);
    F = 1 + 16*(0.53 - E)^3;

    t = mod(4.32e4 * lambda_i + tGPS.seconds, 86400);
    PER = max(sum(beta .* phi_m.^(0:3)), 72000);
    x = 2*pi*(t - 50400)/PER;

    if abs(x) < 1.57
        AMP = max(sum(alpha .* phi_m.^(0:3)), 0);
        delTauG = F * (5e-9 + AMP * (1 - x^2/2 + x^4/24));
    else
        delTauG = F * 5e-9;
    end
    
    % Multiply by scaling factor to correct for frequencies other than L1
    delTauG = (154*10.23e6/fc)^2 * delTauG;
else
    error('Invalid ionosphere model');
end

end