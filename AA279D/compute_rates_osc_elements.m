% rates function for numerical integration for propagating osc elements
% from J2 effects

function [oe_dot] = compute_rates_osc_elements(t, oe)

    mu = 3.986004418e14; % (m^3/s^2)
    J2 = 1.082626e-3;
    Re = 6378137;

    a = oe(1);       % Semi-major axis (m)
    e = oe(2);       % Eccentricity
    i = oe(3);  % Inclination (rad)
    RAAN = oe(4); % Right Ascension of Ascending Node (rad)
    omega = oe(5); % Argument of Perigee (rad)
    M = oe(6); % Mean Anomaly (rad)
    
    e_x = e * cos(omega);
    e_y = e * sin(omega);

    n = sqrt(mu / (a^3)); 

    du_dt = (3/4)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*(sqrt(1 - (e_x^2 + e_y^2))*(3*cos(i)^2) + (5*cos(i)^2 - 1));
    dex_dt = -(3/4)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*e_y*(5*cos(i)^2 - 1);
    dey_dt = (3/4)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*e_x*(5*cos(i)^2 - 1);
    dRAAN_dt = -(3/2)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*cos(i);
    
    oe_dot = zeros(4,1);
    oe_dot(1) = du_dt;
    oe_dot(2) = dex_dt;
    oe_dot(3) = dey_dt;
    oe_dot(4) = dRAAN_dt;
end

