clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

% initial cheif elements
a_TSX_init    = 6886536.686;  % m
e_TSX_init    = 0.0001264;
i_TSX_init    = deg2rad(97.4453);
RAAN_TSX_init = deg2rad(351.0108);
omega_TSX_init= deg2rad(101.2452);
M_TSX_init    = deg2rad(11.6520);
nu_TSX_init   = mean2true(M_TSX_init, e_TSX_init, tol);
u_TSX_init    = M_TSX_init + omega_TSX_init;

scenario = input( ...
    ['Select reconfiguration scenario:\n' ...
     '  1: Test Reconfiguration\n' ...
     '  2: M-D2 to M-D3\n' ...
     '  3: M-D3 to M-D4\n' ...
     'Your choice: '] );

switch scenario
    case 1
        rel_qns_pre = [0, 0, 200, 200, 200, 200];
        rel_qns_post  = [0, 0, 100, 100, 100, 100];
        scenario_name = 'Test Recongifuration';
        
    case 2
        rel_qns_pre = [0, 0, 0, 300, 0, 400];
        rel_qns_post  = [0, 0, 0, 300, 0, 500];
        scenario_name = 'M-D2 to M-D3';

    case 3
        % M-D3 to M-D4
        rel_qns_pre = [0, 0, 0, 300, 0, 500];
        rel_qns_post  = [0, 0, 0, 500, 0, 300];
        scenario_name = 'M-D3 to M-D4';
        
    otherwise
        error('Invalid scenario');
end

fprintf('Running scenario "%s": rel_qns_init = [%g %g %g %g %g %g]\n rel_qns_fin = [%g %g %g %g %g %g]\n', ...
        scenario_name, [rel_qns_pre, rel_qns_post]);

delta_nom = [rel_qns_post(1), rel_qns_post(3:6), 0]./a_TSX_init;
da_nom = rel_qns_post(1);
dlambda_nom  = rel_qns_post(2);

TSX_init_oe = [a_TSX_init, e_TSX_init, i_TSX_init, ...
               RAAN_TSX_init, omega_TSX_init, M_TSX_init];
TSX_init_rv = oe2rv(TSX_init_oe, mu);

TDX_init_oe = qns2oe(TSX_init_oe, rel_qns_pre);
TDX_init_rv = oe2rv(TDX_init_oe, mu);
[TDX_init_rtn, ~] = eci2rtn(TSX_init_rv, TDX_init_rv);
rel_state_init = [TDX_init_rtn; TSX_init_rv];

% timing parameters
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init^3);
T          = 2*pi/n;
n_orbit    = 15;
tend       = n_orbit*T;
num_points = 10000;
dt         = (tend - tstart)/(num_points-1);
t_grid     = linspace(tstart, tend, num_points).';
t_orbit    = t_grid / T;
orbit_num = floor(t_orbit) + 1;

%Lyapunov parameters and thrust limit
k     = 1e3;      % Lyapunov scaling
N_ip  = 14;       % in-plane exponent
N_oop = 14;       % out-of-plane exponent
u_max = 1e-4;     % maximum thrust accel (m/s^2)

%propagation with delta Vs
state_out = zeros(num_points, 12);
state_out(1,:) = rel_state_init';

TSX_oe = zeros(num_points, 6);
TSX_oe(1,:) = TSX_init_oe;
TDX_oe = zeros(num_points, 6);
TDX_oe(1,:) = TDX_init_oe;

rel_oe = zeros(num_points, 6);
rel_oe(1,:) = TSX_init_oe(1)*compute_roes(TSX_init_oe, TDX_init_oe);

TSX_ECI_hist = zeros(num_points,3);
TDX_ECI_hist = zeros(num_points,3);

state_cur = rel_state_init;
TSX_oe_cur = TSX_init_oe;
TDX_oe_cur = TDX_init_oe;
TSX_ECI_cur = TSX_init_rv;
TDX_ECI_cur = TDX_init_rv;
TDX_RTN_cur = TDX_init_rtn;
rel_oe_cur = rel_qns_pre;

TSX_ECI_hist(1,:) = TSX_init_rv(1:3)';    % chief initial ECI pos
TDX_ECI_hist(1,:) = TDX_init_rv(1:3)';  

%Control history
u_hist   = zeros(3, num_points);
phi_hist = zeros(num_points, 1);
a_hist   = zeros(num_points,1);
dv_hist  = zeros(3, num_points);

for idx = 2:num_points

    t_cur = t_grid(idx-1);
    t_next = t_grid(idx);

    %delta V computations
    a_cur = TSX_oe_cur(1);
    n_cur = sqrt(mu/a_cur^3);
    delta_cur = [rel_oe_cur(1), rel_oe_cur(3:6), 0]./a_cur;

    [A_c, B_c] = plant_reduced_qns(TSX_oe_cur);
    A5 = A_c(1:5,1:5);
    B5 = B_c(1:5,:);

    Delta5 = delta_cur(1:5) - delta_nom(1:5);

    phi_ip   = atan2(delta_cur(3), delta_cur(2));
    phi_oop  = atan2(delta_cur(5), delta_cur(4));
    phi_loc  = wrapTo2Pi(TSX_oe_cur(5) + mean2true(TSX_oe_cur(6), TSX_oe_cur(2), tol));
    Jp = phi_loc - phi_ip;
    Hp = phi_loc - phi_oop;
    P5 = (1/k) * diag([cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Hp)^N_oop; cos(Hp)^N_oop]);

    %dv for Delta da for dlambda_dot
    e_cur = TSX_oe_cur(2);
    f_cur = mean2true(TSX_oe_cur(6), e_cur, tol);
    da_cur = rel_oe_cur(1);
    dlambda_cur = rel_oe_cur(2);
    dlambda_dot_cur = -3/2*n_cur*da_cur;
    % 
    Ddlambda = dlambda_cur - dlambda_nom;
    tau = 1000000;
    dlambda_dot_des = - Ddlambda/tau;
    % 
    Dda_tan = -2/3* dlambda_dot_des / n_cur;
    % 
    dvt = Dda_tan * n_cur / (2 * (1+e_cur*cos(f_cur)));
    u_dvt = dvt / dt;

    % control
    u = [0.0; -pinv(B5) * (A5*delta_cur(1:5)' + P5*Delta5')]; % no r dv
    u(2) = u(2) + u_dvt;
    dv = u .* dt;
    un = norm(u);

    % record
    u_hist(:,idx)  = u;
    phi_hist(idx)  = phi_loc;
    a_hist(idx)    = un;
    dv_hist(:, idx)   = dv;

    %apply delta V to RTN
    TDX_RTN_cur(4:6) = TDX_RTN_cur(4:6) + dv;

    % write the modified RTN back into state and fill histories
    [TDX_ECI_cur, ~] = rtn2eci(TSX_ECI_cur, TDX_RTN_cur);

    [t_out, TDX_ECI_next] = ode4(@compute_rates_rv_perturbed, [t_cur, t_next]', TDX_ECI_cur,  dt);
    TDX_ECI_next = TDX_ECI_next(2,:)';
    [t_out, TSX_ECI_next] = ode4(@compute_rates_rv_perturbed, [t_cur, t_next]', TSX_ECI_cur,  dt);
    TSX_ECI_next = TSX_ECI_next(2,:)';

    TSX_ECI_hist(idx,:) = TSX_ECI_next(1:3)'; 
    TDX_ECI_hist(idx,:) = TDX_ECI_next(1:3)';

    [TDX_RTN_next, ~] = eci2rtn(TSX_ECI_next, TDX_ECI_next);
    state_next = [TDX_RTN_next; TSX_ECI_next]';
    state_out(idx,:) = state_next;

    TSX_params = rv2oe(TSX_ECI_next, mu);
    TSX_oe_next = [TSX_params(1:5), true2mean(TSX_params(6), TSX_params(2))];
    TSX_oe_next(3:6) = wrapTo2Pi(TSX_oe_next(3:6));
    TSX_oe(idx,:) = TSX_oe_next;

    TDX_params = rv2oe(TDX_ECI_next, mu);
    TDX_oe_next = [TDX_params(1:5), true2mean(TDX_params(6), TDX_params(2))];
    TDX_oe_next(3:6) = wrapTo2Pi(TDX_oe_next(3:6));
    TDX_oe(idx,:) = TDX_oe_next;

    rel_oe_next = TSX_oe_next(1)*compute_roes(TSX_oe_next, TDX_oe_next)';
    % if (abs(rel_oe_next(2))>10^5)
    %     rel_oe_next(2) = 0.0;
    % end
    rel_oe(idx,:) = rel_oe_next;

    %update current states
    state_cur = state_next;
    TSX_oe_cur = TSX_oe_next;
    TDX_oe_cur = TDX_oe_next;
    TSX_ECI_cur = TSX_ECI_next;
    TDX_ECI_cur = TDX_ECI_next;
    TDX_RTN_cur = TDX_RTN_next;
    rel_oe_cur = rel_oe_next;

end

%Plotting
% --- RTN frame plots ---
% TR projection
figure;
subplot(1,3,1);
plot(state_out(:,2), state_out(:,1));
xlabel('T [m]');
ylabel('R [m]');
grid on;
axis equal;

% NR projection
subplot(1,3,2);
plot(state_out(:,3), state_out(:,1));
xlabel('N [m]');
ylabel('R [m]');
grid on;
axis equal;

% TN projection
subplot(1,3,3);
plot(state_out(:,2), state_out(:,3));
xlabel('T [m]');
ylabel('N [m]');
grid on;
axis equal;

% 3D RTN trajectory
figure;
plot3(state_out(:,1), state_out(:,2), state_out(:,3));
xlabel('R [m]');
ylabel('T [m]');
zlabel('N [m]');
grid on;
axis equal;

% --- Relative orbital elements vs. orbit number ---
figure;
% Δa
subplot(3,2,1);
plot(t_orbit, rel_oe(:,1), 'LineWidth', 1.5);
xlabel('Orbit Number');
ylabel('a\deltaa [m]');
grid on;
% Δ\lambda
subplot(3,2,2);
plot(t_orbit, rel_oe(:,2), 'LineWidth', 1.5);
xlabel('Orbit Number');
ylabel('a\delta\lambda [m]');
grid on;
% Δe_x
subplot(3,2,3);
plot(t_orbit, rel_oe(:,3), 'LineWidth', 1.5);
xlabel('Orbit Number');
ylabel('a\deltae_x [m]');
grid on;
% Δe_y
subplot(3,2,4);
plot(t_orbit, rel_oe(:,4), 'LineWidth', 1.5);
xlabel('Orbit Number');
ylabel('a\deltae_y [m]');
grid on;
% Δi_x
subplot(3,2,5);
plot(t_orbit, rel_oe(:,5), 'LineWidth', 1.5);
xlabel('Orbit Number');
ylabel('a\deltai_x [m]');
grid on;
% Δi_y
subplot(3,2,6);
plot(t_orbit, rel_oe(:,6), 'LineWidth', 1.5);
xlabel('Orbit Number');
ylabel('a\deltai_y [m]');
grid on;

% --- 2D relative‐element scatter for control visualization ---
figure;
% Δλ vs Δa
subplot(1,3,1);
plot(rel_oe(:,2), rel_oe(:,1), 'LineWidth',1.2);
xlabel('a\delta\lambda [m]');
ylabel('a\deltaa [m]');
grid on;
axis equal;
% Δe_x vs Δe_y
subplot(1,3,2);
plot(rel_oe(:,3), rel_oe(:,4), 'LineWidth',1.2);
xlabel('a\deltae_x [m]');
ylabel('a\deltae_y [m]');
grid on;
axis equal;
% Δi_x vs Δi_y
subplot(1,3,3);
plot(rel_oe(:,5), rel_oe(:,6), 'LineWidth',1.2);
xlabel('a\deltai_x [m]');
ylabel('a\deltai_y [m]');
grid on;
axis equal;

% 8) Plot control acceleration level
figure;
plot(t_orbit, a_hist);
xlabel('Orbit Number'); ylabel('Acceleration magnitude (m/s^2)');
%title('Control Acceleration Level');

% assume phi_hist and a_hist are both (N×1)
phi = wrapTo2Pi(phi_hist);
a   = a_hist;

% find where φ jumps backwards by more than π
dphi = diff(phi);
wrapIdx = find(dphi < -pi);

% break the line at each wrap
for k = 1:numel(wrapIdx)
  i = wrapIdx(k)+1;      % index of the wrapped point
  phi(i) = NaN;
  a(i)   = NaN;
end

% now plot without the wrap‐around line
figure
plot(phi, a);
xlabel('Argument of Latitude \phi (rad)')
ylabel('Acceleration magnitude (m/s^2)')
%title('Control vs Argument of Latitude (no wrap)')
grid on

phi_unw = unwrap(phi_hist);
figure
plot(phi_unw, a_hist)
xlabel('Unwrapped Argument of Latitude \phi (rad)')
ylabel('Acceleration magnitude (m/s^2)')
%title('Control vs Unwrapped Argument of Latitude')
grid on

% cumulative delta-v
cum_dv_r = cumsum(abs(dv_hist(1,:)));
cum_dv_t = cumsum(abs(dv_hist(2,:)));
cum_dv_n = cumsum(abs(dv_hist(3,:)));
cum_dv = cum_dv_r + cum_dv_t + cum_dv_n;
total_dv = sum(abs(dv_hist(:)));  
fprintf('Total Δv = %.3f m/s\n', total_dv);
figure;
total_dv = sum(cum_dv);
plot(t_orbit, cum_dv);
xlabel('Orbit Number'); ylabel('Cumulative Delta-v (m/s)');
%title('Cumulative Delta-v');

%plotting absolute orbital elements
figure;
subplot(3,2,1)
plot(t_orbit, TSX_oe(:,1)); grid on;
xlabel('Orbits'); ylabel('a [m]')
title('Semi–major axis')

subplot(3,2,2)
plot(t_orbit, TSX_oe(:,2)); grid on;
xlabel('Orbits'); ylabel('e')
title('Eccentricity')

subplot(3,2,3)
plot(t_orbit, rad2deg(TSX_oe(:,3))); grid on;
xlabel('Orbits'); ylabel('i [deg]')
title('Inclination')

subplot(3,2,4)
plot(t_orbit, rad2deg(TSX_oe(:,4))); grid on;
xlabel('Orbits'); ylabel('\Omega [deg]')
title('RAAN')

subplot(3,2,5)
plot(t_orbit, rad2deg(TSX_oe(:,5))); grid on;
xlabel('Orbits'); ylabel('\omega [deg]')
title('Argument of perigee')

subplot(3,2,6)
plot(t_orbit, rad2deg(TSX_oe(:,6))); grid on;
xlabel('Orbits'); ylabel('M [deg]')
title('Mean anomaly')

figure;
subplot(3,2,1)
plot(t_orbit, TDX_oe(:,1)); grid on;
xlabel('Orbits'); ylabel('a [m]')
title('Semi–major axis')

subplot(3,2,2)
plot(t_orbit, TDX_oe(:,2)); grid on;
xlabel('Orbits'); ylabel('e')
title('Eccentricity')

subplot(3,2,3)
plot(t_orbit, rad2deg(TDX_oe(:,3))); grid on;
xlabel('Orbits'); ylabel('i [deg]')
title('Inclination')

subplot(3,2,4)
plot(t_orbit, rad2deg(TDX_oe(:,4))); grid on;
xlabel('Orbits'); ylabel('\Omega [deg]')
title('RAAN')

subplot(3,2,5)
plot(t_orbit, rad2deg(TDX_oe(:,5))); grid on;
xlabel('Orbits'); ylabel('\omega [deg]')
title('Argument of perigee')

subplot(3,2,6)
plot(t_orbit, rad2deg(TDX_oe(:,6))); grid on;
xlabel('Orbits'); ylabel('M [deg]')
title('Mean anomaly')
title("TSX Orbital Elements");