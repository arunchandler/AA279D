function [rv_dot] = compute_rates_rv_perturbed(t, rv)

    mu = 3.986004418e14; % (m^3/s^2)
    J2 = 1.082626e-3;
    Re = 6378137; % m

    rv_dot = zeros(6, 1);
    r = rv(1:3);  % Position vector (m)
    v = rv(4:6);  % Velocity vector (m/s)
    r_mag = norm(r);
    
    % Compute the acceleration due to gravity
    a_gravity = -mu .* r / r_mag^3;
    
    % Compute J2 perturbation acceleration
    z = r(3);
    coeff_J2 = (3/2) * (mu * J2 * Re^2/r_mag^5);
    ax_J2 = coeff_J2 * (5 * z^2/r_mag^2 - 1) * r(1);
    ay_J2 = coeff_J2 * (5 * z^2/r_mag^2 - 1) * r(2);
    az_J2 = coeff_J2 * (5 * z^2/r_mag^2 - 3) * r(3);

    a_J2 = [ax_J2; ay_J2; az_J2];
    
    % Total acceleration
    a_total = a_gravity + a_J2;

    rv_dot(1:3) = v;
    rv_dot(4:6) = a_total;
end
