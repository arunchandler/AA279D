% rates function for numerical integration for relative motion
% unperturbed
% Inputs:
%   
% Outputs: 

function [rv_dot] = compute_rates_rv_relative_unperturbed(t, )

    mu = 3.986004418e14; % (m^3/s^2)

    r = rv(1:3);
    r_mag = norm(r);
    v = rv(4:6);
    
    rv_dot = zeros(6,1);
    rv_dot(1:3) = v;

    rv_dot(4:6) = -mu .* r / r_mag^3;

end