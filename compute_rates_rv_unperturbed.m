% rates function for numerical integration for propagating rv state
% unperturbed

function [rv_dot] = compute_rates_rv_unperturbed_ECI(t, rv)

    global mu;
    
    rv_dot = zeros(6,1);
    rv_dot(1:3) = rv(4:6);
    rv_dot(4:6) = -mu .* rv(1:3) / norm(rv(1:3))^3;

end