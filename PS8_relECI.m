clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d sigma_pos sigma_vel
tol = 1e-10;
Re  = 6378137;          % m
J2  = 1.082626e-3;
mu  = 3.986004418e14;   % m^3/s^2
s_d = 86400;            % s/day
sigma_pos = 1;          % [m]
sigma_vel = 0.001;       % [m/s]

% initial cheif elements & state
a_TSX_init    = 6886536.686;  % m
e_TSX_init    = 0.0001264;
i_TSX_init    = deg2rad(97.4453);
RAAN_TSX_init = deg2rad(351.0108);
omega_TSX_init= deg2rad(101.2452);
M_TSX_init    = deg2rad(11.6520);
nu_TSX_init   = mean2true(M_TSX_init, e_TSX_init, tol);
u_TSX_init    = M_TSX_init + omega_TSX_init;

TSX_oe_init = [a_TSX_init, e_TSX_init, i_TSX_init, ...
               RAAN_TSX_init, omega_TSX_init, M_TSX_init];
TSX_rv_init = oe2rv(TSX_oe_init, mu);

% initial ROEs, elements, and state for deputy
rel_qns  = [0, 0, 0, 300, 0, 500];
TDX_oe_init = qns2oe(TSX_oe_init, rel_qns);
TDX_rv_init = oe2rv(TDX_oe_init, mu);

rel_rv_init = TDX_rv_init - TSX_rv_init;

% timing parameters
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init^3);
T          = 2*pi/n;
n_orbit    = 15;
tend       = n_orbit*T;
num_points = 1000;
dt         = (tend - tstart)/(num_points-1);
t_grid     = linspace(tstart, tend, num_points).';
t_orbit    = t_grid / T;
orbit_num = floor(t_orbit) + 1;

% Generate ground truth
[t_out, TSX_rv_gt] = ode4(@compute_rates_rv_perturbed, [tstart,tend]', TSX_rv_init, dt);
[t_out, TDX_rv_gt] = ode4(@compute_rates_rv_perturbed, [tstart,tend]', TDX_rv_init, dt);
rel_rv_gt = TDX_rv_gt - TSX_rv_gt;

% Generate measurements (add noise to ground truth)
sigma_rv = diag([sigma_pos^2*ones(3,1);
                  sigma_vel^2*ones(3,1)]);

N = length(t_grid);

noise_TSX = (sqrtm(sigma_rv) * randn(6, N))';
noise_TDX = (sqrtm(sigma_rv) * randn(6, N))';
noise_rel = noise_TDX - noise_TSX;

TSX_rv_meas = TSX_rv_gt + noise_TSX;
TDX_rv_meas = TDX_rv_gt + noise_TDX;
rel_rv_meas = TDX_rv_meas - TSX_rv_meas;

% Set initial estimate and covariance
TSX_x0 = TSX_rv_init + (sqrtm(sigma_rv) * randn(6, 1));
x0 = rel_rv_init + (sqrtm(sigma_rv) * randn(6, 1));
P0 = sigma_rv;

% Define process and measurements noise covariances
Q = 2 * P0;
R = 0.5 * P0; % Might need to change

% Set up and propagate EKF
TSX_rv_ekf = zeros(num_points, 6);
TSX_rv_ekf(1,:) = TSX_x0;
rel_rv_ekf = zeros(num_points, 6);
rel_rv_ekf(1,:) = x0;
TSX_x_k1k1 = TSX_x0;
x_k1k1 = x0;
P_TSX_hist = zeros(num_points,6,6);
P_TSX_hist(1,:,:) = P0;
P_TSX_k1k1 = P0;
P_hist = zeros(num_points,6,6);
P_hist(1,:,:) = P0;
P_k1k1 = P0;
y_TSX_prefit   = zeros(num_points,6);
y_TSX_postfit  = zeros(num_points,6);
y_prefit   = zeros(num_points,6);
y_postfit  = zeros(num_points,6);

for idx = 2:num_points

    t = t_grid(idx-1);
    t_next = t+dt;

    % -------- PREDICT ---------
    % Uncomment for STM propagation ---------
    F = compute_F_rel_j2(x_k1k1(1:3), TSX_x_k1k1(1:3));
    Phi = expm(F*dt); % For STM Propagation
    x_kk1 = Phi * x_k1k1 + (sqrtm(sigma_rv) * randn(6, 1));
    % ---------------------------------------

    % Uncomment for ode propagation ----------
    [~, TSX_x_kk1] = ode4(@compute_rates_rv_perturbed, [t,t_next]', TSX_x_k1k1, dt);
    TSX_x_kk1 = TSX_x_kk1(2,:)' + (sqrtm(sigma_rv) * randn(6, 1));
    % ----------------------------------------

    F = compute_F_rel_j2(x_kk1(1:3), TSX_x_kk1(1:3)); % This F is for continuous EKF, so next lines will be Pdot, not P
    F_TSX = compute_F_j2(TSX_x_kk1);
    Pdot_kk1 = F * P_k1k1 * F.' + Q;
    Pdot_TSX_kk1 = F_TSX * P_k1k1 * F_TSX.' + Q;
    P_kk1 = P_k1k1 + Pdot_kk1 .* dt;
    P_TSX_kk1 = P_TSX_k1k1 + Pdot_TSX_kk1 .* dt;

    y = rel_rv_meas(idx,:)' - x_kk1;
    y_TSX = TSX_rv_meas(idx,:)' - TSX_x_kk1;
    y_prefit(idx,:) = y';
    y_TSX_prefit(idx,:) = y_TSX';

    % -------- UPDATE --------
    y = rel_rv_meas(idx,:)' - x_kk1;
    y_TSX = TSX_rv_meas(idx,:)' - TSX_x_kk1;

    H = eye(6);
    S = H * P_kk1 * H.' + R;
    S_TSX = H * P_TSX_kk1 * H.' + R;

    K = P_kk1 * H.' * inv(S);
    K_TSX = P_TSX_kk1 * H.' * inv(S_TSX);
    
    x_kk = x_kk1 + K * y;
    TSX_x_kk = TSX_x_kk1 + K_TSX * y_TSX;

    P_kk = (eye(6) - K * H) * P_kk1;
    P_TSX_kk = (eye(6) - K_TSX * H) * P_TSX_kk1;
    P_hist(idx, :,:) = P_kk;
    P_TSX_hist(idx, :,:) = P_TSX_kk;

    y_postfit(idx,:) = (rel_rv_meas(idx,:)' - x_kk)';
    y_TSX_postfit(idx,:) = (TSX_rv_meas(idx,:)' - TSX_x_kk)';

    % Store state and update P and x
    rel_rv_ekf(idx,:) = x_kk;
    TSX_rv_ekf(idx,:) = TSX_x_kk;
    x_k1k1 = x_kk;
    TSX_x_k1k1 = TSX_x_kk;
    P_k1k1 = P_kk;
    P_TSX_k1k1 = P_TSX_kk;

end

% Plotting ------------------------------------------
state_labels = {'r_x error [m]','r_y error [m]','r_z error [m]','v_x error [m/s]','v_y error [m/s]','v_z error [m/s]'};

% --- Rel true vs. filtered ---
figure;
for idx = 1:6
    ax_tsx(idx) = subplot(3,2,idx);
    plot(t_orbit, rel_rv_gt(:, idx),  'b-'); hold on;
    plot(t_orbit, rel_rv_ekf(:, idx), '--r');
    ylabel(state_labels{idx});
    if idx >= 5, xlabel('Time (s)'); end
    hold off;
end
legend(ax_tsx(1), {'True','Filtered'}, 'Location','best');
%sgtitle('True vs. Filtered State Trajectories');

% --- TSX true vs. filtered ---
figure;
for idx = 1:6
    ax_tsx(idx) = subplot(3,2,idx);
    plot(t_orbit, TSX_rv_gt(:, idx),  'b-'); hold on;
    plot(t_orbit, TSX_rv_ekf(:, idx), '--r');
    ylabel(state_labels{idx});
    if idx >= 5, xlabel('Time (s)'); end
    hold off;
end
legend(ax_tsx(1), {'True','Filtered'}, 'Location','best');
%sgtitle('TSX: True vs. Filtered State Trajectories');

% Relative filtered error - including covariance bounds
figure;
for idx = 1:6
    subplot(3,2,idx)
    hold on;
    
    err    = rel_rv_ekf(:,idx) - rel_rv_gt(:,idx);
    sigma1 = squeeze(sqrt( P_hist(:,idx,idx) ));
    sigma2 = 2*sigma1;
    
    h_err = plot(t_orbit, err);               % error
    h1    = plot(t_orbit, +sigma1, 'r--');          % +1σ
    plot(   t_orbit, -sigma1, 'r--');
    h2    = plot(t_orbit, +sigma2,'m--');           % +2σ
    plot(   t_orbit, -sigma2,'m--');              % –2σ
    
    ylabel([state_labels{idx} ' error']);
    if idx>4, xlabel('Orbit #'); end
    
    if idx==1
      legend([h_err h1 h2], ...
             'error','\pm1\sigma','\pm2\sigma', ...
             'Location','best');
    end
    
    hold off;
end

% TSX Filtered error - including covariance bounds
figure;
for idx = 1:6
    subplot(3,2,idx)
    hold on;
    
    err    = TSX_rv_ekf(:,idx) - TSX_rv_gt(:,idx);
    sigma1 = squeeze(sqrt( P_TSX_hist(:,idx,idx) ));
    sigma2 = 2*sigma1;
    
    h_err = plot(t_orbit, err);               % error
    h1    = plot(t_orbit, +sigma1, 'r--');          % +1σ
    plot(   t_orbit, -sigma1, 'r--');
    h2    = plot(t_orbit, +sigma2,'m--');           % +2σ
    plot(   t_orbit, -sigma2,'m--');              % –2σ
    
    ylabel([state_labels{idx} ' error']);
    if idx>4, xlabel('Orbit #'); end
    
    if idx==1
      legend([h_err h1 h2], ...
             'error','\pm1\sigma','\pm2\sigma', ...
             'Location','best');
    end
    
    hold off;
end

% True statistics
last_orbit_idx = t_orbit >= (n_orbit-1);  

err_rel = rel_rv_ekf   - rel_rv_gt;    % [num_points×6]
err_TSX = TSX_rv_ekf - TSX_rv_gt;

mean_rel = mean( err_rel(last_orbit_idx,:), 1 );    % 1×6
mean_TSX = mean( err_TSX(last_orbit_idx,:), 1 );
std_rel  =   std( err_rel(last_orbit_idx,:),  0,1 ); % 1×6
std_TSX  =   std( err_TSX(last_orbit_idx,:),  0,1 );

fprintf('\nrel steady‐state error over last orbit:\n');
for i=1:6
    fprintf('  %s:   mean = %+8.3e   std = %8.3e\n', ...
            state_labels{i}, mean_rel(i), std_rel(i));
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
    plot(t_orbit, y_prefit(:,i),  'b-'); hold on;
    plot(t_orbit, noise_rel(:,i),     'r--');
    plot(t_orbit, y_postfit(:,i), 'g-');
    ylabel(state_labels{i});
    if i>4, xlabel('Orbit #'); end
    hold off;
end
% one legend on the first subplot only:
legend(ax_rel(1), {'prefit','injected noise','postfit'}, 'Location','best');

figure;
for i = 1:6
    ax_tsx(i) = subplot(3,2,i);
    plot(t_orbit, y_TSX_prefit(:,i),  'b-'); hold on;
    plot(t_orbit, noise_TSX(:,i),     'r--');
    plot(t_orbit, y_TSX_postfit(:,i), 'g-');
    ylabel(state_labels{i});
    if i>4, xlabel('Orbit #'); end
    hold off;
end
% one legend on the first subplot only:
legend(ax_tsx(1), {'prefit','injected noise','postfit'}, 'Location','best');