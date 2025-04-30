function Ks = getYAconstants(rv_rel_RTN, a, e, f0)
% Compute the integration constants for the Yamanaka–Ankersen
%   Tschauner–Hempel solution using the Yamanaka–Ankersen transition matrix.
%
% Inputs:
%   rv_rel_RTN : 6×1 relative state vector in RTN at initial true anomaly f0
%                [x; y; z; vx; vy; vz], where positions are in meters, 
%                velocities are in m/s
%   a          : chief orbit semi-major axis [m]\%   e          : chief orbit eccentricity
%   f0         : initial true anomaly [rad]
%
% Outputs:
%   Ks : 6×1 vector of integration constants for the YA solution

    global mu;

   
    eta  = sqrt(1 - e^2);
    n    = sqrt(mu/a^3);
    k    = 1 + e*cos(f0);
    kp   = -e*sin(f0);
    tau  = 0;  
    % (since we're evaluating at f = f0 for the initial inversion)

    psi_x_1 = (1/k) + (3/2)*(kp*tau);
    psi_x_2 = sin(f0);
    psi_x_3 = cos(f0);
    
    psi_y_1 = (-3/2)*k*tau;
    psi_y_2 = (1 + (1/k))*cos(f0);
    psi_y_3 = -(1 + (1/k))*sin(f0);
    psi_y_4 = (1/k);

    psi_z_5 = (1/k)*sin(f0);
    psi_z_6 = (1/k)*cos(f0);

    psi_x_dot_1 = (kp/2) - (3/2)*(k^2)*(k-1)*tau;
    psi_x_dot_2 = (k^2)*cos(f0);
    psi_x_dot_3 = -(k^2)*sin(f0);

    psi_y_dot_1 = -(3/2)*(k+(k^2)*kp*tau);
    psi_y_dot_2 = -(k^2 + 1)*sin(f0);
    psi_y_dot_3 = -e - (k^2 + 1)*cos(f0);
    psi_y_dot_4 = -kp;

    psi_z_dot_5 = e + cos(f0);
    psi_z_dot_6 = -sin(f0);

    A = [a*(eta^2)*eye(3), zeros(3); zeros(3), (a*n/eta)*eye(3)];
    B = [psi_x_1 psi_x_2 psi_x_3 0 0 0; ...
        psi_y_1 psi_y_2 psi_y_3 psi_y_4 0 0; ...
        0 0 0 0 psi_z_5 psi_z_6; ...
        psi_x_dot_1 psi_x_dot_2 psi_x_dot_3 0 0 0; ...
        psi_y_dot_1 psi_y_dot_2 psi_y_dot_3 psi_y_dot_4 0 0; ...
        0 0 0 0 psi_z_dot_5 psi_z_dot_6];
    total = A * B;
    Ks = total\rv_rel_RTN;
end
