% rates function for numerical integration for propagating osc elements
% from J2 effects

function [oe_dot] = compute_rates_osc_elements(t, oe)

    mu = 3.986004418e14; % (m^3/s^2)
    J2 = 1.082626e-3;
    Re = 6378137;

    a = oe(1);      % Semi-major axis (m)
    e_x = oe(2);     % Eccentricity x component
    e_y = oe(3);     % Eccentricity y component
    i = oe(4);      % Inclination (rad)
    RAAN = oe(5);   % (rad)
    u = oe(6);      % Argument of latitude (rad)

    n = sqrt(mu / (a^3)); 

    da_dt = 0;
    dex_dt = -(3/4)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*e_y*(5*cos(i)^2 - 1);
    dey_dt = (3/4)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*e_x*(5*cos(i)^2 - 1);
    di_dt = 0;
    dRAAN_dt = -(3/2)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*cos(i);
    du_dt = (3/4)*n*J2*(Re/(a*(1 - (e_x^2 + e_y^2))))^2*(sqrt(1 - (e_x^2 + e_y^2))*(3*cos(i)^2) + (5*cos(i)^2 - 1));
    
    oe_dot = zeros(1,6);
    oe_dot(1) = da_dt;
    oe_dot(2) = dex_dt;
    oe_dot(3) = dey_dt;
    oe_dot(4) = di_dt;
    oe_dot(5) = dRAAN_dt;
    oe_dot(6) = du_dt;
end

