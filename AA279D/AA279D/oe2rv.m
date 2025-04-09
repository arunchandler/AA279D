% Computes position and velocity in ECI coordinates from orbital elements
%   Inputs:
%       koe - orbital elements in the form [a,e,i,RAAN,w,M]
%       mu - gravitational parameter in m^3/s^2

function rv_eci = oe2rv(oe, mu)

    a = oe(1);       % Semi-major axis (m)
    e = oe(2);       % Eccentricity
    I = oe(3);  % Inclination (rad)
    RAAN = oe(4); % Right Ascension of Ascending Node (rad)
    omega = oe(5); % Argument of Perigee (rad)
    M = oe(6); % Mean Anomaly (rad)
    
    nu = mean2true(M,e,1e-10);
    
    % Compute Perifocal Coordinates (r, v)
    p = a * (1 - e^2); 
    r_perifocal = [p * cos(nu) / (1 + e * cos(nu));
                   p * sin(nu) / (1 + e * cos(nu));
                   0];
               
    v_perifocal = sqrt(mu / p) * [-sin(nu);
                                   e + cos(nu);
                                   0];

    % Rotation Matrices
    R3_W = [ cos(-RAAN), sin(-RAAN), 0;
            -sin(-RAAN), cos(-RAAN), 0;
             0, 0, 1];

    R1_i = [1, 0, 0;
            0, cos(-I), sin(-I);
            0, -sin(-I), cos(-I)];
        
    R3_w = [ cos(-omega), sin(-omega), 0;
            -sin(-omega), cos(-omega), 0;
             0, 0, 1];

    % Total Rotation Matrix: Perifocal to ECI
    Q_p2eci = R3_W * R1_i * R3_w;

    % Comte ECI Position and Velocity
    r_eci = Q_p2eci * r_perifocal;
    v_eci = Q_p2eci * v_perifocal;
    rv_eci = [r_eci;v_eci];
end
