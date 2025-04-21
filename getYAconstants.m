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

    %-- Orbit parameters --
    p  = a * (1 - e^2);          % semi-latus rectum [m]
    h  = sqrt(mu * p);           % specific angular momentum [m^2/s]

    %-- Chief radius and anomaly-rate at f0 --
    r0  = p / (1 + e * cos(f0)); % radius at f0 [m]
    fdot = h / r0^2;             % df/dt [rad/s]

    %-- Normalize the RTN state --
    pos = rv_rel_RTN(1:3);       % [m]
    vel = rv_rel_RTN(4:6);       % [m/s]

    % normalized position
    x0 = pos(1) / r0;
    y0 = pos(2) / r0;
    z0 = pos(3) / r0;

    % convert velocity to dρ/df and normalize
    drdf = vel ./ fdot;          % derivative w.r.t. true anomaly [m/rad]
    xp0  = drdf(1) / r0;
    yp0  = drdf(2) / r0;
    zp0  = drdf(3) / r0;

    X0 = [x0; y0; z0; xp0; yp0; zp0];

    %-- Build Yamanaka-Ankersen transition matrix at f0 --

    %--- orbit and anomaly parameters at f0 ---
    eta  = sqrt(1 - e^2);
    n    = sqrt(mu/a^3);
    k    = 1 + e*cos(f0);
    kp   = -e*sin(f0);
    tau  = 0;  
    % (since we're evaluating at f = f0 for the initial inversion)

    %--- first- and second-derivatives ---
    k_sin    = k*sin(f0);
    k_cos    = k*cos(f0);
    dk_sin   = cos(f0) + e*cos(2*f0);
    dk_cos   = - ( sin(f0) + e*sin(2*f0) );

    Phi = zeros(6,6);

    % row 1
    Phi(1,1) = 1 + 1.5 * k * kp * tau;
    Phi(1,2) = k_sin;
    Phi(1,3) = k_cos;

    % row 2
    Phi(2,1) = -1.5 * k^2 * tau;
    Phi(2,2) = (1 + k)*cos(f0);
    Phi(2,3) = -(1 + k)*sin(f0);
    Phi(2,4) = 1;

    % row 3
    Phi(3,5) = sin(f0);
    Phi(3,6) = cos(f0);

    % row 4
    Phi(4,1) = 1.5*(kp/k) - 1.5*e*dk_sin*tau;
    Phi(4,2) = dk_sin;
    Phi(4,3) = dk_cos;

    % row 5
    Phi(5,1) = -1.5*(1 + 2*k*kp*tau);
    Phi(5,2) = -2 * k * sin(f0);
    Phi(5,3) = e - 2*k * cos(f0);

    % row 6
    Phi(6,5) = cos(f0);
    Phi(6,6) = -sin(f0);

    %-- invert to get your six K constants --
    Ks = Phi \ X0;
end
