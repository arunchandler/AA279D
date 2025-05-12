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
fprintf('Running "%s" with nominal ROE = [%g %g %g %g %g]\n', ...
        scenario_name, rel_nom);

% append zero for delta_a_dot
delta_nom = [rel_nom; 0];

% 2) Chief initial orbital elements (mean): [a; e; i; RAAN; omega; M]
a0     = 6886536.686;
e0     = 0.0001264;
i0     = deg2rad(97.4453);
RAAN0  = deg2rad(351.0108);
omega0 = deg2rad(101.2452);
M0     = deg2rad(11.6520);
oe_c   = [a0; e0; i0; RAAN0; omega0; M0];

% 3) Time grid
n0       = sqrt(mu/a0^3);
Torbit   = 2*pi/n0;
n_orbit  = 15;
t_end    = n_orbit * Torbit;
num_pts  = 1000;
t_grid   = linspace(0, t_end, num_pts)';
dt       = t_grid(2) - t_grid(1);

% 4) Lyapunov parameters and thrust limit
k     = 1e3;      % Lyapunov scaling
N_ip  = 14;       % in-plane exponent
N_oop = 14;       % out-of-plane exponent
u_max = 1e-4;     % maximum thrust accel (m/s^2)

% 5) Integrate reduced-state ODE
delta0 = delta_nom;   % start at nominal
delta0 = [0; 0; 400; 0; 500; 0]; %for reconfiguration
odefun = @(t,delta) reduced_dynamics(delta, oe_c, delta_nom, k, N_ip, N_oop, u_max);
[~, hist_delta] = ode4(odefun, [0 t_end], delta0, dt);

% 6) Plot state trajectories
orbit_number = t_grid / Torbit;
figure;
plot(orbit_number, hist_delta(:,1:5));
xlabel('Orbit Number'); ylabel('ROE Components');
legend('delta_a','delta_ex','delta_ey','delta_ix','delta_iy','Location','best');
title('Reduced ROE vs Orbits');

% 7) Compute and plot delta-v history
u_hist  = zeros(2, num_pts);
dv_hist = zeros(num_pts,1);
for idx = 1:num_pts
    delta_i = hist_delta(idx,:)';
    [A_c, B_c] = plant_reduced_qns(oe_c);
    Delta     = delta_i - delta_nom;
    phi_ip    = atan2(delta_i(4), delta_i(3));
    phi_oop   = atan2(delta_i(6), delta_i(5));
    phi = oe_c(5) + mean2true(oe_c(6), oe_c(2));
    Jp = phi - phi_ip;
    Hp = phi - phi_oop;
    P         = (1/k)*diag([cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Hp)^N_oop; cos(Hp)^N_oop; 1]);
    u_i       = - pinv(B_c)*(A_c*delta_i + P*Delta);
    u_i       = max(min(u_i, u_max), -u_max);
    u_hist(:,idx) = u_i;
    dv_hist(idx)   = norm(u_i)*dt;
end
cum_dv = cumsum(dv_hist);

figure;
plot(orbit_number, dv_hist);
xlabel('Orbit Number'); ylabel('Delta-v per step (m/s)');
title('Delta-v per Step');

figure;
plot(orbit_number, cum_dv);
xlabel('Orbit Number'); ylabel('Cumulative Delta-v (m/s)');
title('Cumulative Delta-v');

% 8) Plot ROE error norm
figure;
element_labels = {'a\delta a [m]', 'a\deltae_x [m]', 'a\deltae_y [m]', 'a\deltai_x [m]', 'a\deltai_y [m]'};

for i = 1:5
    subplot(3,2,i);
    plot(orbit_number, hist_delta(:,i) - rel_nom(i));
    xlabel('Orbit Number');
    ylabel(element_labels{i});
    title(['ROE Error: ' element_labels{i}]);
    grid on;
end


%% standalone ODE function
function delta_dot = reduced_dynamics(delta, oe_c, delta_nom, k, N_ip, N_oop, u_max)
    %delta = [delta_a; delta_ex; delta_ey; delta_ix; delta_iy; delta_a_dot];
    % 1) reduced plant & input
    [A_c, B_c] = plant_reduced_qns(oe_c);
    % 2) tracking error
    Delta = delta - delta_nom;
    % 3) Lyapunov P
    phi_ip  = atan2(delta(4), delta(3));
    phi_oop = atan2(delta(6), delta(5));
    phi = oe_c(5) + mean2true(oe_c(6), oe_c(2));
    Jp = phi - phi_ip;
    Hp = phi - phi_oop;
    P       = (1/k)*diag([cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Hp)^N_oop; cos(Hp)^N_oop; 1]);
    % 4) control
    u = - pinv(B_c)*(A_c*delta + P*Delta);
    % 5) saturation
    u = max(min(u, u_max), -u_max);
    % 6) ODE
    delta_dot = A_c*delta + B_c*u;
end
