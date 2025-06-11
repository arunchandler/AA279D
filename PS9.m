clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d sigma_pos sigma_vel sigma_qns_a sigma_qns_lam sigma_qns_e sigma_qns_i
tol = 1e-10;
Re  = 6378137;          % m
J2  = 1.082626e-3;
mu  = 3.986004418e14;   % m^3/s^2
s_d = 86400;            % s/day
sigma_pos = 1;
sigma_vel = 0.005;
sigma_qns_a = 3; 
sigma_qns_lam = 10;
sigma_qns_e = 5; 
sigma_qns_i = 5; 

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

fprintf('Running scenario "%s" (with noise): rel_qns_init = [%g %g %g %g %g %g]\n rel_qns_fin = [%g %g %g %g %g %g]\n', ...
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
N = length(t_grid);

%Lyapunov parameters and thrust limit
k     = 1e3;      % Lyapunov scaling
N_ip  = 14;       % in-plane exponent
N_oop = 14;       % out-of-plane exponent
u_max = 1e-4;     % maximum thrust accel (m/s^2)

% ------------------------
% Ground truth propagation
% ------------------------

state_out = zeros(num_points, 12);
state_out(1,:) = rel_state_init';

TSX_oe = zeros(num_points, 6);
TSX_oe(1,:) = TSX_init_oe;
TDX_oe = zeros(num_points, 6);
TDX_oe(1,:) = TDX_init_oe;

rel_oe = zeros(num_points, 6);
rel_oe(1,:) = TSX_init_oe(1)*compute_roes(TSX_init_oe, TDX_init_oe);

TSX_ECI_hist = zeros(num_points,6);
TDX_ECI_hist = zeros(num_points,6);

state_cur = rel_state_init;
TSX_oe_cur = TSX_init_oe;
TDX_oe_cur = TDX_init_oe;
TSX_ECI_cur = TSX_init_rv;
TDX_ECI_cur = TDX_init_rv;
TDX_RTN_cur = TDX_init_rtn;
rel_oe_cur = rel_qns_pre;

TSX_ECI_hist(1,:) = TSX_init_rv';    % chief initial ECI pos
TDX_ECI_hist(1,:) = TDX_init_rv';  

%Control history
u_hist   = zeros(3, num_points);
phi_hist = zeros(num_points, 1);
a_hist   = zeros(num_points,1);
dv_hist  = zeros(3, num_points);

sigma_rv = diag([sigma_pos^2*ones(3,1);
                  sigma_vel^2*ones(3,1)]);

sigma_roe_init = diag([sigma_qns_a^2, sigma_qns_lam^2, sigma_qns_e^2, sigma_qns_e^2, sigma_qns_i^2, sigma_qns_i^2]);

sigma_roe_meas = diag([0.5^2 5^2 1^2 1^2 0.5^2 0.5^2]);

noise_roe = (sqrtm(sigma_roe_meas) * randn(6, N))';
noise_TSX = (sqrtm(sigma_rv) * randn(6, N))';
noise_TDX = (sqrtm(sigma_rv) * randn(6, N))';

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
    % TDX_ECI_next = TDX_ECI_next(2,:)' + noise_TDX(idx,:)';
    TDX_ECI_next = TDX_ECI_next(2,:)';

    [t_out, TSX_ECI_next] = ode4(@compute_rates_rv_perturbed, [t_cur, t_next]', TSX_ECI_cur,  dt);
    % TSX_ECI_next = TSX_ECI_next(2,:)' + noise_TSX(idx,:)';
    TSX_ECI_next = TSX_ECI_next(2,:)';

    TSX_ECI_hist(idx,:) = TSX_ECI_next'; 
    TDX_ECI_hist(idx,:) = TDX_ECI_next';

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

    % rel_oe_next = TSX_oe_next(1)*compute_roes(TSX_oe_next, TDX_oe_next)' + noise_roe(idx,:);
    rel_oe_next = TSX_oe_next(1)*compute_roes(TSX_oe_next, TDX_oe_next)';
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

% cumulative delta-v
dv_mag  = sqrt( dv_hist(1,:).^2 + dv_hist(2,:).^2 + dv_hist(3,:).^2 );
cum_dv   = cumsum( dv_mag );
total_dv = cum_dv(end);
fprintf('Total Δv = %.3f m/s\n', total_dv);
figure;
plot(t_orbit, cum_dv);
xlabel('Orbit Number'); ylabel('Cumulative Delta-v (m/s)');
%title('Cumulative Delta-v');

% --------------------------
% Filtered state propagation
% --------------------------

%Lyapunov parameters and thrust limit
k     = 1e3;      % Lyapunov scaling
N_ip  = 14;       % in-plane exponent
N_oop = 14;       % out-of-plane exponent
u_max = 1e-4;     % maximum thrust accel (m/s^2)

state_out_filt = zeros(num_points, 12);
state_out_filt(1,:) = rel_state_init';

TSX_oe_filt = zeros(num_points, 6);
TSX_oe_filt(1,:) = TSX_init_oe;
TDX_oe_filt = zeros(num_points, 6);
TDX_oe_filt(1,:) = TDX_init_oe;

rel_oe_filt = zeros(num_points, 6);
rel_oe_filt(1,:) = TSX_init_oe(1)*compute_roes(TSX_init_oe, TDX_init_oe);

TSX_ECI_hist_filt = zeros(num_points,6);
TDX_ECI_hist_filt = zeros(num_points,6);

state_cur_filt = rel_state_init;
TSX_oe_cur_filt = TSX_init_oe;
TDX_oe_cur_filt = TDX_init_oe;
TSX_ECI_cur_filt = TSX_init_rv;
TDX_ECI_cur_filt = TDX_init_rv;
TDX_RTN_cur_filt = TDX_init_rtn;
rel_oe_cur_filt = rel_qns_pre;

TSX_ECI_hist_filt(1,:) = TSX_init_rv';
TDX_ECI_hist_filt(1,:) = TDX_init_rv';  

TSX_x0 = TSX_init_rv + (sqrtm(sigma_rv) * randn(6, 1));
x0 = rel_qns_pre' + (sqrtm(sigma_roe_meas) * randn(6, 1));
x0_unfiltered = rel_qns_pre';
P0 = sigma_rv;
P0_roe = sigma_roe_init;

% base Q on the diff between STM and GVE
Q = P0/10;
% measurement covariances
R      = diag([2^2*ones(3,1); (0.02)^2*ones(3,1)]);  % absolute GPS

% process noise (ROE)
n   = sqrt(mu/a_TSX_init^3);
sigma_a_t = 2.5e-4;                      % 40 µm/s²  ← tuned
Q_roe = 100*(sigma_a_t^2 / n^2) * diag([4 4 2 2 1 1]) * dt;
R_roe = P0_roe;

TSX_x_k1k1 = TSX_x0;
x_k1k1 = x0;
P_TSX_hist = zeros(num_points,6,6);
P_TSX_hist(1,:,:) = P0;
P_TSX_k1k1 = P0;
P_hist = zeros(num_points,6,6);
P_hist(1,:,:) = P0_roe;
P_k1k1 = P0_roe;
resid_TSX_prefit   = zeros(num_points,6);
resid_TSX_postfit  = zeros(num_points,6);
resid_prefit   = zeros(num_points,6);
resid_postfit  = zeros(num_points,6);
NIS = zeros(num_points, 1);
NEES = zeros(num_points, 1); 

%Control history
u_hist_filt   = zeros(3, num_points);
phi_hist_filt = zeros(num_points, 1);
a_hist_filt   = zeros(num_points,1);
dv_hist_filt  = zeros(3, num_points);

for idx = 2:num_points

    t_cur_filt = t_grid(idx-1);
    t_next_filt = t_grid(idx);

    %delta V computations
    a_cur_filt = TSX_oe_cur_filt(1);
    n_cur_filt = sqrt(mu/a_cur_filt^3);
    delta_cur_filt = [rel_oe_cur_filt(1), rel_oe_cur_filt(3:6), 0]./a_cur_filt;

    [A_c_filt, B_c_filt] = plant_reduced_qns(TSX_oe_cur_filt);
    A5_filt = A_c_filt(1:5,1:5);
    B5_filt = B_c_filt(1:5,:);

    Delta5_filt = delta_cur_filt(1:5) - delta_nom(1:5);

    phi_ip_filt   = atan2(delta_cur_filt(3), delta_cur_filt(2));
    phi_oop_filt  = atan2(delta_cur_filt(5), delta_cur_filt(4));
    phi_loc_filt  = wrapTo2Pi(TSX_oe_cur_filt(5) + mean2true(TSX_oe_cur_filt(6), TSX_oe_cur_filt(2), tol));
    Jp_filt = phi_loc_filt - phi_ip_filt;
    Hp_filt = phi_loc_filt - phi_oop_filt;
    P5_filt = (1/k) * diag([cos(Jp_filt)^N_ip; cos(Jp_filt)^N_ip; cos(Jp_filt)^N_ip; cos(Hp_filt)^N_oop; cos(Hp_filt)^N_oop]);

    %dv for Delta da for dlambda_dot
    e_cur_filt = TSX_oe_cur_filt(2);
    f_cur_filt = mean2true(TSX_oe_cur_filt(6), e_cur_filt, tol);
    da_cur_filt = rel_oe_cur_filt(1);
    dlambda_cur_filt = rel_oe_cur_filt(2);
    dlambda_dot_cur_filt = -3/2*n_cur_filt*da_cur_filt;
    % 
    Ddlambda_filt = dlambda_cur_filt - dlambda_nom;
    tau = 1000000;
    dlambda_dot_des_filt = - Ddlambda_filt/tau;
    % 
    Dda_tan_filt = -2/3* dlambda_dot_des_filt / n_cur_filt;
    % 
    dvt_filt = Dda_tan_filt * n_cur_filt / (2 * (1+e_cur_filt*cos(f_cur_filt)));
    u_dvt_filt = dvt_filt / dt;

    % control
    u_filt = [0.0; -pinv(B5_filt) * (A5_filt*delta_cur_filt(1:5)' + P5_filt*Delta5_filt')]; % no r dv
    u_filt(2) = u_filt(2) + u_dvt_filt;
    dv_filt = u_filt .* dt;
    un_filt = norm(u_filt);

    % record
    u_hist_filt(:,idx)  = u_filt;
    phi_hist_filt(idx)  = phi_loc_filt;
    a_hist_filt(idx)    = un_filt;
    dv_hist_filt(:, idx)   = dv_filt;

    % -------- PREDICT ---------
    % ROE STM propagation ---------
    F = stm_qns_j2(dt, TSX_oe_cur_filt);
    B = get_B(TSX_oe_cur_filt);
    G = B .* dt;
    x_kk1 = F * x_k1k1 + G * u_filt;
    % -----------------------------

    % TSX ode propagation ----------
    [~, TSX_x_kk1] = ode4(@compute_rates_rv_perturbed, [t_cur_filt,t_next_filt]', TSX_x_k1k1, dt);
    TSX_x_kk1 = TSX_x_kk1(2,:)';

    F_TSX = compute_F_j2(TSX_x_kk1); % This F is for continuous EKF, so next lines will be Pdot, not P
    Pdot_TSX_kk1 = F_TSX * P_TSX_k1k1 + P_TSX_k1k1 * F_TSX.' + Q;
    P_kk1 = F * P_k1k1 * F.' + Q_roe;
    P_TSX_kk1 = P_TSX_k1k1 + Pdot_TSX_kk1 .* dt;

    resid = rel_oe(idx,:)'+ noise_roe(idx,:)' - x_kk1;
    resid_TSX = TSX_ECI_hist(idx,:)'+noise_TSX(idx,:)' - TSX_x_kk1;
    resid_prefit(idx,:) = resid';
    resid_TSX_prefit(idx,:) = resid_TSX';

    % -------- UPDATE --------
    y = resid;
    y_TSX = resid_TSX;

    H = eye(6);
    S = H * P_kk1 * H.' + R_roe;
    S_TSX = H * P_TSX_kk1 * H.' + R;

    K = P_kk1 * H.' * inv(S);
    K_TSX = P_TSX_kk1 * H.' * inv(S_TSX);

    x_kk = x_kk1 + K * y;
    TSX_x_kk = TSX_x_kk1 + K_TSX * y_TSX;

    P_kk = (eye(6) - K * H) * P_kk1;
    P_TSX_kk = (eye(6) - K_TSX * H) * P_TSX_kk1;
    P_hist(idx, :,:) = P_kk;
    P_TSX_hist(idx, :,:) = P_TSX_kk;

    resid_postfit(idx,:) = rel_oe(idx,:)+noise_roe(idx,:) - x_kk';
    resid_TSX_postfit(idx,:) = TSX_ECI_hist(idx,:)+noise_TSX(idx,:) - TSX_x_kk';

    nu  = resid_prefit(idx,:)';      % innovation r_k
    NIS(idx)  = nu.' / S * nu;       % χ²(6) should lie in [1.6,14.4] 95 %
    e   = x_kk - rel_oe(idx,:)';
    NEES(idx) = e.'  / P_kk * e;   % χ²(6) should lie in [0.9,16.8] 95 %

    % write the modified RTN back into state and fill histories
    TSX_ECI_next_filt = TSX_x_kk;
    rel_oe_next_filt = x_kk';
    TSX_params = rv2oe(TSX_ECI_next_filt, mu);
    TSX_oe_next_filt = [TSX_params(1:5), true2mean(TSX_params(6), TSX_params(2))];
    TSX_oe_next_filt(3:6) = wrapTo2Pi(TSX_oe_next_filt(3:6));
    TDX_oe_next_filt = qns2oe(TSX_oe_next_filt, rel_oe_next_filt);
    TDX_ECI_next_filt = oe2rv(TDX_oe_next_filt, mu);

    TSX_ECI_hist_filt(idx,:) = TSX_ECI_next_filt; 
    TDX_ECI_hist_filt(idx,:) = TDX_ECI_next_filt;

    [TDX_RTN_next_filt, ~] = eci2rtn(TSX_ECI_next_filt, TDX_ECI_next_filt);
    state_next_filt = [TDX_RTN_next_filt; TSX_ECI_next_filt]';
    state_out_filt(idx,:) = state_next_filt;

    TSX_oe_filt(idx,:) = TSX_oe_next_filt;
    TDX_oe_filt(idx,:) = TDX_oe_next_filt;
    rel_oe_filt(idx,:) = rel_oe_next_filt;

    %update current states
    state_cur_filt = state_next_filt;
    TSX_oe_cur_filt = TSX_oe_next_filt;
    TDX_oe_cur_filt = TDX_oe_next_filt;
    TSX_ECI_cur_filt = TSX_ECI_next_filt;
    TDX_ECI_cur_filt = TDX_ECI_next_filt;
    TDX_RTN_cur_filt = TDX_RTN_next_filt;
    rel_oe_cur_filt = rel_oe_next_filt;

    x_k1k1 = x_kk;
    TSX_x_k1k1 = TSX_x_kk;
    P_k1k1 = P_kk;
    P_TSX_k1k1 = P_TSX_kk;

end

% Plotting ------------------------------------------
state_labels_error = {'r_x error [m]','r_y error [m]','r_z error [m]','v_x error [m/s]','v_y error [m/s]','v_z error [m/s]'};
roe_labels_error = {'a\deltaa error [m]','a\delta\lambda error [m]','a\deltae_x error [m]','a\deltae_y error [m/s]','a\deltai_x error [m/s]','a\deltai_y error [m/s]'};

roe_labels = {'a\deltaa [m]','a\delta\lambda [m]','a\deltae_x [m]','a\deltae_y [m]','a\deltai_x [m]','a\deltai_y [m]'};
state_labels = {'r_x [m]','r_y [m]','r_z [m]','v_x [m/s]','v_y [m/s]','v_z [m/s]'};

% --- Rel true vs. filtered ---
figure;
for idx = 1:6
    ax_rel(idx) = subplot(3,2,idx);
    hold on;

    % time, truth, estimate, and 1σ
    t     = t_orbit;
    truth = rel_oe(:,idx);
    est   = rel_oe_filt(:,idx);
    sigma = sqrt( squeeze( P_hist(:,idx,idx) ) );

    % build the closed‐loop for shading (est+σ then est–σ)
    x2       = [t; flipud(t)];
    band1    = [ est + sigma; flipud(est - sigma) ];
    fill(x2, band1, [1 0.6 0.6], ...
         'EdgeColor','none','FaceAlpha',0.4);

    % now plot truth and filtered mean on top
    plot(t, est,   '--r');
    plot(t, truth, 'b-');

    ylabel(roe_labels{idx});
    if idx >= 5, xlabel('Orbit #'); end

    hold off;
end

% legend only needs to go on the first subplot
legend(ax_rel(1), {'±1σ bound','Filtered', 'True'}, 'Location','best');
%legend(ax_rel(1), {'True','Filtered','Unfiltered'}, 'Location','best');
%sgtitle('True vs. Filtered State Trajectories');

% --- TSX true vs. filtered ---
figure;
for idx = 1:6
    ax_tsx(idx) = subplot(3,2,idx);
    plot(t_orbit, TSX_ECI_hist(:, idx),  'b-'); hold on;
    plot(t_orbit, TSX_ECI_hist_filt(:, idx), '--r');
    ylabel(state_labels{idx});
    if idx >= 5, xlabel('Orbit #'); end
    hold off;
end
legend(ax_tsx(1), {'True','Filtered'}, 'Location','best');
%sgtitle('TSX: True vs. Filtered State Trajectories');

% Relative filtered error - including covariance bounds
figure;
for idx = 1:6
    subplot(3,2,idx)
    hold on;
    
    err    = rel_oe_filt(:,idx) - rel_oe(:,idx);
    sigma1 = squeeze(sqrt( P_hist(:,idx,idx) ));
    sigma3 = 3*sigma1;
    
    h_err = plot(t_orbit, err);               % error
    h1    = plot(t_orbit, +sigma1, 'r--');          % +1σ
    plot(   t_orbit, -sigma1, 'r--');
    h2    = plot(t_orbit, +sigma3,'m--');           % +2σ
    plot(   t_orbit, -sigma3,'m--');              % –2σ
    
    ylabel([roe_labels_error{idx}]);
    if idx>4, xlabel('Orbit #'); end
    
    if idx==1
      legend([h_err h1 h2], ...
             'error','\pm1\sigma','\pm3\sigma', ...
             'Location','best');
    end
    
    hold off;
end

% TSX Filtered error - including covariance bounds
figure;
for idx = 1:6
    subplot(3,2,idx)
    hold on;
    
    err    = TSX_ECI_hist_filt(:,idx) - TSX_ECI_hist(:,idx);
    sigma1 = squeeze(sqrt( P_TSX_hist(:,idx,idx) ));
    sigma3 = 3*sigma1;
    
    h_err = plot(t_orbit, err);               % error
    h1    = plot(t_orbit, +sigma1, 'r--');          % +1σ
    plot(   t_orbit, -sigma1, 'r--');
    h2    = plot(t_orbit, +sigma3,'m--');           % +3σ
    plot(   t_orbit, -sigma3,'m--');              % –3σ
    
    ylabel([state_labels_error{idx}]);
    if idx>4, xlabel('Orbit #'); end
    
    if idx==1
      legend([h_err h1 h2], ...
             'error','\pm1\sigma','\pm3\sigma', ...
             'Location','best');
    end
    
    hold off;
end

% True statistics
last_orbit_idx = t_orbit >= (n_orbit-1);  

err_rel = rel_oe_filt   - rel_oe;    % [num_points×6]
err_TSX = TSX_ECI_hist_filt - TSX_ECI_hist;

mean_rel = mean( err_rel(last_orbit_idx,:), 1 );    % 1×6
mean_TSX = mean( err_TSX(last_orbit_idx,:), 1 );
std_rel  =   std( err_rel(last_orbit_idx,:),  0,1 ); % 1×6
std_TSX  =   std( err_TSX(last_orbit_idx,:),  0,1 );

fprintf('\nRel steady‐state error over last orbit:\n');
for i=1:6
    fprintf('  %s:   mean = %+8.3e   std = %8.3e\n', ...
            roe_labels{i}, mean_rel(i), std_rel(i));
end

fprintf('\nTSX steady‐state error over last orbit:\n');
for i=1:6
    fprintf('  %s:   mean = %+8.3e   std = %8.3e\n', ...
            state_labels{i}, mean_TSX(i), std_TSX(i));
end

% Pre-fit & Post-fit residuals
figure;
for i = 1:6
    ax_rel(i) = subplot(3,2,i);
    plot(t_orbit, resid_prefit(:,i),  'b-'); hold on;
    plot(t_orbit, noise_roe(:,i),     'r--');
    plot(t_orbit, resid_postfit(:,i), 'g-');
    ylabel(roe_labels{i});
    if i>4, xlabel('Orbit #'); end
    hold off;
end
% one legend on the first subplot only:
legend(ax_rel(1), {'prefit','injected noise','postfit'}, 'Location','best');

figure;
for i = 1:6
    ax_tsx(i) = subplot(3,2,i);
    plot(t_orbit, resid_TSX_prefit(:,i),  'b-'); hold on;
    plot(t_orbit, noise_TSX(:,i),     'r--');
    plot(t_orbit, resid_TSX_postfit(:,i), 'g-');
    ylabel(state_labels{i});
    if i>4, xlabel('Orbit #'); end
    hold off;
end
% one legend on the first subplot only:
legend(ax_tsx(1), {'prefit','injected noise','postfit'}, 'Location','best');

avg_NEES = mean(NEES);

% Print it
fprintf('\nAverage NEES over all steps: %8.3e\n', avg_NEES);

%Plotting
% --- RTN frame plots ---
% TR projection
figure;
subplot(1,3,1);
plot(state_out_filt(:,2), state_out_filt(:,1));
xlabel('T [m]');
ylabel('R [m]');
grid on;
axis equal;

% NR projection
subplot(1,3,2);
plot(state_out_filt(:,3), state_out_filt(:,1));
xlabel('N [m]');
ylabel('R [m]');
grid on;
axis equal;

% TN projection
subplot(1,3,3);
plot(state_out_filt(:,2), state_out_filt(:,3));
xlabel('T [m]');
ylabel('N [m]');
grid on;
axis equal;

% 3D RTN trajectory
figure;
plot3(state_out_filt(:,1), state_out_filt(:,2), state_out_filt(:,3));
xlabel('R [m]');
ylabel('T [m]');
zlabel('N [m]');
grid on;
axis equal;

% --- 2D relative‐element scatter for control visualization ---
figure;
% Δλ vs Δa
subplot(1,3,1);
plot(rel_oe_filt(:,2), rel_oe_filt(:,1), 'LineWidth',1.2);
xlabel('a\delta\lambda [m]');
ylabel('a\deltaa [m]');
grid on;
axis equal;
% Δe_x vs Δe_y
subplot(1,3,2);
plot(rel_oe_filt(:,3), rel_oe_filt(:,4), 'LineWidth',1.2);
xlabel('a\deltae_x [m]');
ylabel('a\deltae_y [m]');
grid on;
axis equal;
% Δi_x vs Δi_y
subplot(1,3,3);
plot(rel_oe_filt(:,5), rel_oe_filt(:,6), 'LineWidth',1.2);
xlabel('a\deltai_x [m]');
ylabel('a\deltai_y [m]');
grid on;
axis equal;

% 8) Plot control acceleration level
figure;
plot(t_orbit, a_hist_filt);
xlabel('Orbit Number'); ylabel('Acceleration magnitude (m/s^2)');
%title('Control Acceleration Level');

% assume phi_hist and a_hist are both (N×1)
phi_filt = wrapTo2Pi(phi_hist_filt);
a_filt   = a_hist_filt;

% find where φ jumps backwards by more than π
dphi_filt = diff(phi_filt);
wrapIdx = find(dphi_filt < -pi);

% break the line at each wrap
for k = 1:numel(wrapIdx)
  i = wrapIdx(k)+1;      % index of the wrapped point
  phi_filt(i) = NaN;
  a_filt(i)   = NaN;
end

% now plot without the wrap‐around line
figure
plot(phi_filt, a);
xlabel('Argument of Latitude \phi (rad)')
ylabel('Acceleration magnitude (m/s^2)')
%title('Control vs Argument of Latitude (no wrap)')
grid on

% cumulative delta-v
dv_mag_filt  = sqrt( dv_hist_filt(1,:).^2 + dv_hist_filt(2,:).^2 + dv_hist_filt(3,:).^2 );
cum_dv_filt   = cumsum( dv_mag_filt );
total_dv_filt = cum_dv_filt(end);
fprintf('Total Δv = %.3f m/s\n', total_dv_filt);
figure;
total_dv = sum(cum_dv_filt);
plot(t_orbit, cum_dv_filt);
xlabel('Orbit Number'); ylabel('Cumulative Delta-v (m/s)');
%title('Cumulative Delta-v');