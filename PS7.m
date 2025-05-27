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

% Generate ground truth
[t_out, TSX_rv_gt] = ode4(@compute_rates_rv_perturbed, [tstart,tend]', TSX_rv_init, dt);
[t_out, TDX_rv_gt] = ode4(@compute_rates_rv_perturbed, [tstart,tend]', TDX_rv_init, dt);

% Generate measurements (add noise to ground truth)
sigma_rv = diag([sigma_pos^2*ones(3,1);
                  sigma_vel^2*ones(3,1)]);

N = length(t_grid);

noise_TSX = (sqrtm(sigma_rv) * randn(6, N))';
noise_TDX = (sqrtm(sigma_rv) * randn(6, N))';

TSX_rv_meas = TSX_rv_gt + noise_TSX;
TDX_rv_meas = TDX_rv_gt + noise_TDX;

% Set initial estimate and covariance
TSX_x0 = TSX_rv_init + (sqrtm(sigma_rv) * randn(6, 1));
TDX_x0 = TDX_rv_init + (sqrtm(sigma_rv) * randn(6, 1));
P0 = sigma_rv;

% Define process and measurements noise covariances
Q = 2 * P0;
R = 0.5 * P0; % Might need to change

% Set up and propagate EKF
TSX_rv_ekf_prefit = zeros(num_points,6);
TDX_rv_ekf_prefit = zeros(num_points,6);
TSX_rv_ekf_prefit(1,:) = TSX_x0;
TDX_rv_ekf_prefit(1,:) = TDX_x0;
TSX_rv_ekf = zeros(num_points, 6);
TDX_rv_ekf = zeros(num_points, 6);
TSX_rv_ekf(1,:) = TSX_x0;
TDX_rv_ekf(1,:) = TDX_x0;
TSX_x_k1k1 = TSX_x0;
TDX_x_k1k1 = TDX_x0;
P_TSX_k1k1 = P0;
P_TDX_k1k1 = P0;

for idx = 2:num_points

    t = t_grid(idx-1);
    t_next = t+dt;

    % -------- UPDATE ---------
    % Uncomment for STM propagation ---------
    % F_TSX = compute_F_j2(TSX_x_k1k1);
    % F_TDX = compute_F_j2(TDX_x_k1k1);
    % Phi_TSX = expm(F_TSX*dt); % For STM Propagation
    % Phi_TDX = expm(F_TDX*dt);
    % TSX_x_kk1 = Phi_TSX * TSX_x_k1k1 + (sqrtm(sigma_rv) * randn(6, 1));
    % TDX_x_kk1 = Phi_TDX * TDX_x_k1k1 + (sqrtm(sigma_rv) * randn(6, 1));
    % ---------------------------------------

    % Uncomment for ode propagation ----------
    [~, TSX_x_kk1] = ode4(@compute_rates_rv_perturbed, [t,t_next]', TSX_x_k1k1, dt);
    [~, TDX_x_kk1] = ode4(@compute_rates_rv_perturbed, [t,t_next]', TDX_x_k1k1, dt);
    TSX_x_kk1 = TSX_x_kk1(2,:)' + (sqrtm(sigma_rv) * randn(6, 1));
    TDX_x_kk1 = TDX_x_kk1(2,:)' + (sqrtm(sigma_rv) * randn(6, 1));
    % ----------------------------------------

    F_TSX = compute_F_j2(TSX_x_kk1(1:3)); % This F is for continuous EKF, so next lines will be Pdot, not P
    F_TDX = compute_F_j2(TDX_x_kk1(1:3));
    Pdot_TSX_kk1 = F_TSX * P_TSX_k1k1 * F_TSX.' + Q;
    Pdot_TDX_kk1 = F_TDX * P_TDX_k1k1 * F_TDX.' + Q;
    P_TSX_kk1 = P_TSX_k1k1 + Pdot_TSX_kk1 .* dt;
    P_TDX_kk1 = P_TDX_k1k1 + Pdot_TDX_kk1 .* dt;

    TSX_rv_ekf_prefit(idx, :) = TSX_x_kk1;
    TDX_rv_ekf_prefit(idx, :) = TDX_x_kk1;

    % -------- PREDICT --------
    y_TSX = TSX_rv_meas(idx,:)' - TSX_x_kk1;
    y_TDX = TDX_rv_meas(idx,:)' - TDX_x_kk1;

    H_TSX = eye(6);
    H_TDX = eye(6);
    S_TSX = H_TSX * P_TSX_kk1 * H_TSX.' + R;
    S_TDX = H_TDX * P_TDX_kk1 * H_TDX.' + R;

    K_TSX = P_TSX_kk1 * H_TSX.' * inv(S_TSX);
    K_TDX = P_TDX_kk1 * H_TDX.' * inv(S_TDX);
    
    TSX_x_kk = TSX_x_kk1 + K_TSX * y_TSX;
    TDX_x_kk = TDX_x_kk1 + K_TDX * y_TDX;

    P_TSX_kk = (eye(6) - K_TSX * H_TSX) * P_TSX_kk1;
    P_TDX_kk = (eye(6) - K_TDX * H_TDX) * P_TDX_kk1;

    % Store state and update P and x
    TSX_rv_ekf(idx,:) = TSX_x_kk;
    TDX_rv_ekf(idx,:) = TDX_x_kk;
    TSX_x_k1k1 = TSX_x_kk;
    TDX_x_k1k1 = TSX_x_kk;
    P_TSX_k1k1 = P_TSX_kk;
    P_TDX_k1k1 = P_TDX_kk;

end

% Plotting ------------------------------------------
state_labels = {'r_x error [m]','r_y error [m]','r_z error [m]','v_x error [m/s]','v_y error [m/s]','v_z error [m/s]'};

%True & filtered states
figure;
for idx = 1:6
    subplot(3,2,idx)
    plot(t_orbit, TSX_rv_gt(:, idx), '-')
    hold on
    plot(t_orbit, TSX_rv_ekf(:, idx), '--')
    ylabel(state_labels{idx})
    if idx == 5 || idx == 6
        xlabel('Time (s)')
    end
end

figure;
for idx = 1:6
    subplot(3,2,idx)
    plot(t_orbit, TDX_rv_gt(:, idx), '-')
    hold on
    plot(t_orbit, TDX_rv_ekf(:, idx), '--')
    ylabel(state_labels{idx})
    if idx == 5 || idx == 6
        xlabel('Time (s)')
    end
end

figure;
for idx = 1:6
    subplot(3,2,idx)
    plot(t_orbit, TSX_rv_ekf_prefit(:, idx) - TSX_rv_gt(:, idx), '-')
    ylabel(state_labels{idx})
    if idx == 5 || idx == 6
        xlabel('Time (s)')
    end
end

figure;
for idx = 1:6
    subplot(3,2,idx)
    plot(t_orbit, TSX_rv_ekf(:, idx) - TSX_rv_gt(:, idx), '-')
    ylabel(state_labels{idx})
    if idx == 5 || idx == 6
        xlabel('Time (s)')
    end
end

% figure;
% plot3(TSX_rv_gt(:,1), TSX_rv_gt(:,2), TSX_rv_gt(:,3), 'k-', 'DisplayName','GT'); 
% hold on;
% plot3(TSX_rv_ekf(:,1), TSX_rv_ekf(:,2), TSX_rv_ekf(:,3), 'b--', 'DisplayName','EKF');
% grid on;    axis equal;
% xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
% legend('Location','best');
% title('TSX: 3D Trajectory—EKF vs. Ground Truth');