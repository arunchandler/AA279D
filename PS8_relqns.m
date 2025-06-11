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
               RAAN_TSX_init, omega_TSX_init, M_TSX_init]';
TSX_rv_init = oe2rv(TSX_oe_init, mu);

% initial ROEs, elements, and state for deputy
rel_qns_init  = [0, 0, 0, 300, 0, 500]';
TDX_oe_init = qns2oe(TSX_oe_init, rel_qns_init);
TDX_rv_init = oe2rv(TDX_oe_init, mu);

% using GVEs cause they maintain numerical stability for longer

% timing parameters
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init^3);
disp(n)
T          = 2*pi/n;
n_orbit    = 5;
tend       = n_orbit*T;
num_points = 1000;
dt         = (tend - tstart)/(num_points-1);
t_grid     = linspace(tstart, tend, num_points).';
t_orbit    = t_grid / T;
orbit_num = floor(t_orbit) + 1;
disp(dt)

% Generate ground truth
[~, TSX_rv_gt] = ode4(@compute_rates_rv_perturbed, [tstart,tend]', TSX_rv_init, dt);
[~, TDX_rv_gt] = ode4(@compute_rates_rv_perturbed, [tstart,tend]', TDX_rv_init, dt);
[~, TSX_oe_gt] = ode4(@compute_rates_GVE_J2, [tstart,tend]', TSX_oe_init, dt);
[~, TDX_oe_gt] = ode4(@compute_rates_GVE_J2, [tstart,tend]', TDX_oe_init, dt);
rel_oe_gt = zeros(num_points,6);
for i = 1:num_points
    rel_oe_gt(i,:) = TSX_oe_gt(i,1)*compute_roes(TSX_oe_gt(i,:), TDX_oe_gt(i,:));
end

% Generate measurements (add noise to ground truth)
sigma_rv = diag([sigma_pos^2*ones(3,1);
                  sigma_vel^2*ones(3,1)]);

sigma_roe_init = diag([sigma_qns_a^2, sigma_qns_lam^2, sigma_qns_e^2, sigma_qns_e^2, sigma_qns_i^2, sigma_qns_i^2]);

sigma_roe_meas = diag([0.5^2 5^2 1^2 1^2 0.5^2 0.5^2]);

N = length(t_grid);

noise_TSX = (sqrtm(sigma_rv) * randn(6, N))';
noise_TDX = (sqrtm(sigma_rv) * randn(6, N))';
noise_roe = (sqrtm(sigma_roe_meas) * randn(6, N))';

TSX_rv_meas = TSX_rv_gt + noise_TSX;
TDX_rv_meas = TDX_rv_gt + noise_TDX;
roe_meas = rel_oe_gt + noise_roe;

% Set initial estimate and covariance
TSX_x0 = TSX_rv_init + (sqrtm(sigma_rv) * randn(6, 1));
x0 = rel_qns_init + (sqrtm(sigma_roe_meas) * randn(6, 1));
x0_unfiltered = rel_qns_init;
P0 = sigma_rv;
P0_roe = sigma_roe_init;

% use differences in RTN as measurements and convert that to ROE
% Error should not leave the covariance bounds, which means we are
% confident

% Define process and measurements noise covariances

% base Q on the diff between STM and GVE
Q = P0/10;
% measurement covariances
R      = diag([2^2*ones(3,1); (0.02)^2*ones(3,1)]);  % absolute GPS

% process noise (ROE)
n   = sqrt(mu/a_TSX_init^3);
sigma_a_t = 2.5e-4;                      % 40 µm/s²  ← tuned
Q_roe = 100*(sigma_a_t^2 / n^2) * diag([4 4 2 2 1 1]) * dt;
R_roe = P0_roe;

disp(Q_roe/R_roe)
% delta lambda hardest to observe but if feeding the whole state there
% shouldn't be drift


% numerical propagation vs STM osculating (PSET 4)
roe_unfiltered = zeros(num_points,6);
roe_unfiltered(1,:) = x0_unfiltered;
for i = 2:num_points
    t = t_grid(i-1);
    Phi = stm_qns_j2(dt, TSX_oe_gt(i-1,:));
    x0_unfiltered = Phi * x0_unfiltered;
    roe_unfiltered(i,:) = x0_unfiltered;
end

%——— Plot relative orbital elements ———%
figure;

% 1) Δa vs Δλ
subplot(3,1,1);
plot(rel_oe_gt(:,1),  rel_oe_gt(:,2),  'b-',  'LineWidth',1.5); hold on;
plot(roe_unfiltered(:,1), roe_unfiltered(:,2), 'r--','LineWidth',1.5);
xlabel('a\deltaa [m]'); ylabel('a\delta\lambda [m]');
axis equal;
legend('Ground truth','Unfiltered STM','Location','best');

% 2) e_x vs e_y
subplot(3,1,2);
plot(rel_oe_gt(:,3),  rel_oe_gt(:,4),  'b-',  'LineWidth',1.5); hold on;
plot(roe_unfiltered(:,3), roe_unfiltered(:,4), 'r--','LineWidth',1.5);
xlabel('a\deltae_x [m]'); ylabel('a\deltae_y [m]');
axis equal;

% 3) i_x vs i_y
subplot(3,1,3);
plot(rel_oe_gt(:,5),  rel_oe_gt(:,6),  'b-',  'LineWidth',1.5); hold on;
plot(roe_unfiltered(:,5), roe_unfiltered(:,6), 'r--','LineWidth',1.5);
xlabel('a\deltai_x [m]'); ylabel('a\deltai_y [m]');
axis equal;

% Set up and propagate EKF
TSX_rv_ekf = zeros(num_points, 6);
TSX_oe_ekf = zeros(num_points, 6);
TSX_oe_ekf(1,:) = TSX_oe_init;
TSX_rv_ekf(1,:) = TSX_x0;
roe_ekf = zeros(num_points, 6);
roe_ekf(1,:) = x0;
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

TSX_oe = TSX_oe_init;

for idx = 2:num_points

    t = t_grid(idx-1);
    t_next = t+dt;

    % -------- PREDICT ---------
    % ROE STM propagation ---------
    F = stm_qns_j2(dt, TSX_oe);
    x_kk1 = F * x_k1k1;
    % -----------------------------

    % TSX ode propagation ----------
    [~, TSX_x_kk1] = ode4(@compute_rates_rv_perturbed, [t,t_next]', TSX_x_k1k1, dt);
    TSX_x_kk1 = TSX_x_kk1(2,:)';
    % ----------------------------------------
    
    F_TSX = compute_F_j2(TSX_x_kk1); % This F is for continuous EKF, so next lines will be Pdot, not P
    Pdot_TSX_kk1 = F_TSX * P_TSX_k1k1 + P_TSX_k1k1 * F_TSX.' + Q;
    P_kk1 = F * P_k1k1 * F.' + Q_roe;
    P_TSX_kk1 = P_TSX_k1k1 + Pdot_TSX_kk1 .* dt;

    resid = roe_meas(idx,:)' - x_kk1;
    resid_TSX = TSX_rv_meas(idx,:)' - TSX_x_kk1;
    resid_prefit(idx,:) = resid';
    resid_TSX_prefit(idx,:) = resid_TSX';

    % -------- UPDATE --------
    y = resid;
    y_TSX = resid_TSX;
    % H should be from ROE to ECI

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

    resid_postfit(idx,:) = (roe_meas(idx,:)' - x_kk)';
    resid_TSX_postfit(idx,:) = (TSX_rv_meas(idx,:)' - TSX_x_kk)';

    % Store state and update P and x

    nu  = resid_prefit(idx,:)';      % innovation r_k
    NIS(idx)  = nu.' / S * nu;       % χ²(6) should lie in [1.6,14.4] 95 %
    e   = x_kk - rel_oe_gt(idx,:)';
    NEES(idx) = e.'  / P_kk * e;   % χ²(6) should lie in [0.9,16.8] 95 %
    % print a header once per iteration (optional)
    %fprintf('--- Step %d ---\n', idx);

    % now print NIS and NEES with the values
    %fprintf('NIS = %+8.3e   NEES = %8.3e\n', ...
           % NIS(idx), NEES(idx));


    roe_ekf(idx,:) = x_kk;
    TSX_rv_ekf(idx,:) = TSX_x_kk;
    x_k1k1 = x_kk;
    TSX_x_k1k1 = TSX_x_kk;
    P_k1k1 = P_kk;
    P_TSX_k1k1 = P_TSX_kk;

    TSX_params = rv2oe(TSX_x_kk, mu);
    TSX_oe = osc2mean([TSX_params(1:5), true2mean(TSX_params(6), TSX_params(2))]);


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
    truth = rel_oe_gt(:,idx);
    est   = roe_ekf(:,idx);
    sigma = sqrt( squeeze( P_hist(:,idx,idx) ) );

    % build the closed‐loop for shading (est+σ then est–σ)
    x2       = [t; flipud(t)];
    band1    = [ est + sigma; flipud(est - sigma) ];
    fill(x2, band1, [1 0.6 0.6], ...
         'EdgeColor','none','FaceAlpha',0.4);

    % now plot truth and filtered mean on top
    plot(t, truth, 'b-', 'LineWidth', 1.5);
    plot(t, est,   '--r', 'LineWidth', 1.5);

    ylabel(roe_labels{idx});
    if idx >= 5, xlabel('Orbit #'); end

    hold off;
end

% legend only needs to go on the first subplot
legend(ax_rel(1), {'±1σ bound','True','Filtered'}, 'Location','best');
%legend(ax_rel(1), {'True','Filtered','Unfiltered'}, 'Location','best');
%sgtitle('True vs. Filtered State Trajectories');

% --- TSX true vs. filtered ---
figure;
for idx = 1:6
    ax_tsx(idx) = subplot(3,2,idx);
    plot(t_orbit, TSX_rv_gt(:, idx),  'b-'); hold on;
    plot(t_orbit, TSX_rv_ekf(:, idx), '--r');
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
    
    err    = roe_ekf(:,idx) - rel_oe_gt(:,idx);
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
    
    err    = TSX_rv_ekf(:,idx) - TSX_rv_gt(:,idx);
    sigma1 = squeeze(sqrt( P_TSX_hist(:,idx,idx) ));
    sigma3 = 3*sigma1;
    
    h_err = plot(t_orbit, err);               % error
    h1    = plot(t_orbit, +sigma1, 'r--');          % +1σ
    plot(   t_orbit, -sigma1, 'r--');
    h2    = plot(t_orbit, +sigma3,'m--');           % +2σ
    plot(   t_orbit, -sigma3,'m--');              % –2σ
    
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

err_rel = roe_ekf   - rel_oe_gt;    % [num_points×6]
err_TSX = TSX_rv_ekf - TSX_rv_gt;

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