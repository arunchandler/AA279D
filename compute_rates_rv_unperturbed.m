% rates function for numerical integration for propagating rv state
% unperturbed

function [rv_dot] = compute_rates_rv_unperturbed(t, rv)

    mu = 3.986004418e14; % (m^3/s^2)
    
    rv_dot = zeros(6,1);
    rv_dot(1:3) = rv(4:6);
    rv_dot(4:6) = -mu .* rv(1:3) / norm(rv(1:3))^3;


end