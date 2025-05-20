%% main_reduced_control.m
% Continuous formation-keeping via reduced-state ROE ODE + Lyapunov feedback
clc; clear; close all;
addpath('mean_osc','control_funcs');   % add your utility function paths
format long g;

global mu J2 Re tol s_d
mu   = 3.986004418e14;    % m^3/s^2
J2   = 1.082626e-3;
Re   = 6378137;           % m
tol  = 1e-10;
s_d  = 86400;             % s/day

% 1) Scenario and nominal reduced ROE (no mean-longitude)
scenario = input(...
  'Select formation-keeping scenario:\n 1: Hold M-D2\n 2: Hold M-D3\nYour choice: ');
switch scenario
  case 1
    rel_nom = [0; 0; 300; 0; 400];  % [delta_a; delta_ex; delta_ey; delta_ix; delta_iy]
    scenario_name = 'Hold M-D2';
  case 2
    rel_nom = [0; 0; 300; 0; 500];
    scenario_name = 'Hold M-D3';
  otherwise
    error('Invalid scenario');
end
fprintf('Running "%s" with nominal ROE = [%g %g %g %g %g]\n', scenario_name, rel_nom);

% append zero for delta_a_dot
delta_nom = [rel_nom; 0];

% 2) Chief initial orbital elements (mean): [a; e; i; RAAN; omega; M]
a0     = 6886536.686;
e0     = 0.0001264;
i0     = deg2rad(97.4453);
RAAN0  = deg2rad(351.0108);
omega0 = deg2rad(101.2452);
M0     = deg2rad(11.6520);
oe_c0  = [a0; e0; i0; RAAN0; omega0; M0];

% 3) Time grid
n0       = sqrt(mu/a0^3);
Torbit   = 2*pi/n0;
n_orbit  = 15;
t_end    = n_orbit * Torbit;
num_pts  = 10000;
t_grid   = linspace(0, t_end, num_pts)';
dt       = t_grid(2) - t_grid(1);

% 4) Lyapunov parameters and thrust limit
k     = 1e3;      % Lyapunov scaling
N_ip  = 14;       % in-plane exponent
N_oop = 14;       % out-of-plane exponent
u_max = 1e-4;     % maximum thrust accel (m/s^2)

% 5) Integrate augmented ODE (ROE-error + OE)
delta0 = [0; 0; 400; 0; 500; 0];
x0 = [delta0; oe_c0];
odefun = @(t,x) augmented_dynamics(t, x, delta_nom, k, N_ip, N_oop, u_max);
[~, hist_x] = ode4(odefun, [0 t_end], x0, dt);

% split histories
hist_delta = hist_x(:,1:6);
hist_oe    = hist_x(:,7:12);

% 6) Plot individual ROE error components vs orbit
orbit_number = t_grid / Torbit;
figure;
components = {'delta_a','delta_ex','delta_ey','delta_ix','delta_iy','delta_a_dot'};
for j = 1:6
    subplot(3,2,j);
    plot(orbit_number, hist_delta(:,j));
    xlabel('Orbit Number'); ylabel(components{j});
    title(components{j});
end
sgtitle('Reduced ROE Error Components vs Orbits');

% 7) Compute control acceleration history & record phi
u_hist   = zeros(2, num_pts);
phi_hist = zeros(num_pts, 1);
a_hist   = zeros(num_pts,1);
dv_hist  = zeros(num_pts,1);
dir_hist = zeros(num_pts,1);
for idx = 1:num_pts
    delta_i = hist_delta(idx,:)';
    oe_i    = hist_oe(idx,:)';

    % plant and gain shaping at this OE
    [A_c, B_c] = plant_reduced_qns(oe_i);
    A5 = A_c(1:5,1:5);
    B5 = B_c(1:5,:);

    % tracking error
    Delta5 = delta_i(1:5) - rel_nom;

    % phases
    phi_ip   = atan2(delta_i(3), delta_i(2));
    phi_oop  = atan2(delta_i(5), delta_i(4));
    phi_loc  = wrapTo2Pi(oe_i(5) + mean2true(oe_i(6), oe_i(2), tol));

    % P5
    Jp = phi_loc - phi_ip;
    Hp = phi_loc - phi_oop;
    P5 = (1/k) * diag([cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Hp)^N_oop; cos(Hp)^N_oop]);

    % control
    u = - pinv(B5) * (A5*delta_i(1:5) + P5*Delta5);

    % record
    u_hist(:,idx)  = u;
    phi_hist(idx)  = phi_loc;
    a_hist(idx)    = norm(u);
    dv_hist(idx)   = a_hist(idx)*dt;
    dir_hist(idx)  = atan2(u(2), u(1));
end
cum_dv = cumsum(dv_hist);

% 8) Plot control acceleration level
figure;
plot(orbit_number, a_hist);
xlabel('Orbit Number'); ylabel('Acceleration magnitude (m/s^2)');
title('Control Acceleration Level');

figure;
plot(phi_hist, u_hist');
xlabel('Argument of Latitude (rad)'); ylabel('Acceleration (m/s^2)');
title('Control vs Argument of Latitude');

% thrust direction
figure;
plot(orbit_number, dir_hist);
xlabel('Orbit Number'); ylabel('Thrust Direction (rad)');
title('Thrust Direction Angle');

% cumulative delta-v
figure;
plot(orbit_number, cum_dv);
xlabel('Orbit Number'); ylabel('Cumulative Delta-v (m/s)');
title('Cumulative Delta-v');

%% augmented_dynamics.m
function x_dot = augmented_dynamics(t, x, delta_nom, k, N_ip, N_oop, u_max)
    global tol;
    % x = [delta(6×1); oe(6×1)]
    delta = x(1:6);
    oe    = x(7:12);

    % 1) propagate OE under J2
    oe_dot = compute_rates_GVE_J2(t, oe);

    % 2) build plant
    [A_c, B_c] = plant_reduced_qns(oe);

    % 3) Lyapunov shaping
    delta5   = delta(1:5);
    Delta5   = delta5 - delta_nom(1:5);
    phi_ip   = atan2(delta5(3), delta5(2));
    phi_oop  = atan2(delta5(5), delta5(4));
    phi_full = oe(5) + mean2true(oe(6), oe(2), tol);

    Jp = phi_full - phi_ip;
    Hp = phi_full - phi_oop;
    P5 = (1/k) * diag([cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Hp)^N_oop; cos(Hp)^N_oop]);

    % 4) compute control
    A5 = A_c(1:5,1:5);
    B5 = B_c(1:5,:);
    u  = - pinv(B5) * (A5*delta5 + P5*Delta5);
    % optional: saturate u here

    % 5) propagate delta under thrust
    delta_dot = A_c*delta + B_c*u;

    % pack
    x_dot = [delta_dot; oe_dot];
end
