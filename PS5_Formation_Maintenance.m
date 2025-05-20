clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
% Constants
tol    = 1e-9;
Re     = 6378137;          % m
J2     = 1.082626e-3;
mu     = 3.986004418e14;   % m^3/s^2
s_d    = 86400;            % seconds per day

%% Scenario selection
scenario = input([
    'Select reconfiguration scenario:\n',...
    '  1: Test Reconfiguration\n',...
    '  2: M-D2 to M-D3\n',...
    'Your choice: ']);
switch scenario
    case 1
        rel_qns_pre  = [0, 0, 200, 200, 200, 200];
        rel_qns_post = [0, 0, 200, 200, 200, 200];
        scenario_name = 'Test Reconfiguration';
    case 2
        rel_qns_pre  = [0, 0,   0, 300,   0, 400];
        rel_qns_post = [0, 0,   0, 300,   0, 400];
        scenario_name = 'M-D2 to M-D3';
    otherwise
        error('Invalid scenario');
end
fprintf('Running scenario "%s"\n', scenario_name);
fprintf('  Pre  ROE: [%g %g %g %g %g %g]\n', rel_qns_pre);
fprintf('  Post ROE: [%g %g %g %g %g %g]\n', rel_qns_post);

%% Chief initial orbital elements and state
a_TSX_init    = 6886536.686;      % m
e_TSX_init    = 0.0001264;
i_TSX_init    = deg2rad(97.4453);
RAAN_TSX_init = deg2rad(351.0108);
omega_TSX_init= deg2rad(101.2452);
M_TSX_init    = deg2rad(11.6520);
TSX_init_oe   = [a_TSX_init, e_TSX_init, i_TSX_init, RAAN_TSX_init, omega_TSX_init, M_TSX_init];
TSX_init_rv   = oe2rv(TSX_init_oe, mu);

%% Deputy initial orbit via QNS
TDX_init_oe = qns2oe(TSX_init_oe, rel_qns_pre);
TDX_init_rv = oe2rv(TDX_init_oe, mu);
[TDX_init_rtn, ~] = eci2rtn(TSX_init_rv, TDX_init_rv);
rel_state_init = [TDX_init_rtn; TSX_init_rv];

%% Time grid & preallocate
n       = sqrt(mu/a_TSX_init^3);
T       = 2*pi/n;
n_orbit = 15;
tstart  = 0;
tend    = n_orbit * T;
num_points = 10000;
dt      = (tend - tstart)/(num_points-1);
t_grid  = linspace(tstart, tend, num_points)';
t_orbit = t_grid / T;
orbit_num = floor(t_orbit) + 1;

state_out     = zeros(num_points, 12);
state_out(1,:) = rel_state_init';

TSX_oe_hist   = zeros(num_points,6);
TSX_oe_hist(1,:) = TSX_init_oe;

TDX_oe_hist   = zeros(num_points,6);
TDX_oe_hist(1,:) = TDX_init_oe;

rel_oe_hist   = zeros(num_points,6);
rel_oe_hist(1,:) = TSX_init_oe(1) * compute_roes(TSX_init_oe, TDX_init_oe)';

TSX_ECI_hist  = zeros(num_points,3);
TDX_ECI_hist  = zeros(num_points,3);

TSX_ECI_cur   = TSX_init_rv;
TDX_ECI_cur   = TDX_init_rv;

TSX_ECI_hist(1,:) = TSX_init_rv(1:3)';
TDX_ECI_hist(1,:) = TDX_init_rv(1:3)';

TDX_RTN_cur   = TDX_init_rtn;

%% Impulsive station-keeping setup
a_chief      = a_TSX_init;
roe_nom      = (rel_qns_post') ./ a_chief;     % desired normalized ROE

delta_de_max = 5 / a_chief;
delta_di_max = 5 / a_chief;

dphi = asin(delta_de_max / norm(roe_nom(3:4)));
de_des = [roe_nom(3)*cos(dphi) - roe_nom(4)*sin(dphi);
          roe_nom(3)*sin(dphi) + roe_nom(4)*cos(dphi)];

dv_times = [];
dv_vals  = [];
dv_hist  = zeros(3, num_points);

%% Main propagation + control loop
for idx = 2:num_points
    t_cur  = t_grid(idx-1);
    t_next = t_grid(idx);
    
    % 1) Compute current normalized ROE
    roe_cur = compute_roes(TSX_oe_hist(idx-1,:), TDX_oe_hist(idx-1,:));
    
    % 2) Detect e-plane drift
    if isempty(dv_vals) && norm(roe_cur(3:4) - roe_nom(3:4)) > delta_de_max
        roe_des       = roe_nom;
        roe_des(3:4) = de_des;
        t_imp = linspace(t_cur, t_cur+T, 5);
        t_imp = t_imp(2:end);
        delta_v_seq = naive_least_squares( ...
            t_imp, roe_cur, roe_des, ...
            TSX_oe_hist(idx-1,:), wrapTo2Pi(TSX_oe_hist(idx-1,5)+TSX_oe_hist(idx-1,6)), ...
            t_cur, t_cur+T);
        dv_times = arrayfun(@(t)find(abs(t_grid-t)==min(abs(t_grid-t)),1), t_imp);
        dv_vals  = delta_v_seq;
    end
    
    % 3) Detect i-plane drift
    if isempty(dv_vals) && norm(roe_cur(5:6) - roe_nom(5:6)) > delta_di_max
        roe_des       = roe_nom;
        roe_des(5:6)= roe_nom(5:6) - delta_di_max*sign(roe_cur(5));
        t_imp = linspace(t_cur, t_cur+T, 3);
        t_imp = t_imp(2:end);
        delta_v_seq = naive_least_squares( ...
            t_imp, roe_cur, roe_des, ...
            TSX_oe_hist(idx-1,:), wrapTo2Pi(TSX_oe_hist(idx-1,5)+TSX_oe_hist(idx-1,6)), ...
            t_cur, t_cur+T);
        dv_times = arrayfun(@(t)find(abs(t_grid-t)==min(abs(t_grid-t)),1), t_imp);
        dv_vals  = delta_v_seq;
    end

    hit = find(dv_times==idx,1);
    if ~isempty(hit)
        dV_rtn = dv_vals(hit,:)';            % [Δv_R; Δv_T; Δv_N] in m/s
        TDX_RTN_cur(4:6) = TDX_RTN_cur(4:6) + dV_rtn;  % bump the RTN-velocities
        % convert back to ECI position+velocity
        [TDX_ECI_cur, ~] = rtn2eci(TSX_ECI_cur, TDX_RTN_cur);
        dv_hist(:,idx) = dV_rtn;
        dv_vals(hit,:)  = [];
        dv_times(hit)   = [];
    end
    
    % 5) Propagate both satellites one step via ODE4
    [~,TDX_out] = ode4(@compute_rates_rv_perturbed, [t_cur t_next]', TDX_ECI_cur, dt);
    [~,TSX_out] = ode4(@compute_rates_rv_perturbed, [t_cur t_next]', TSX_ECI_cur, dt);
    TDX_ECI_next = TDX_out(end,:)';
    TSX_ECI_next = TSX_out(end,:)';
    
    % 6) Update RTN, OE, histories
    [TDX_RTN_next,~] = eci2rtn(TSX_ECI_next, TDX_ECI_next);
    TSX_params = rv2oe(TSX_ECI_next, mu);
    TDX_params = rv2oe(TDX_ECI_next, mu);
    TSX_oe_next = [TSX_params(1:5), true2mean(TSX_params(6),TSX_params(2))];
    TDX_oe_next = [TDX_params(1:5), true2mean(TDX_params(6),TDX_params(2))];
    
    TSX_ECI_hist(idx,:) = TSX_ECI_next(1:3)';
    TDX_ECI_hist(idx,:) = TDX_ECI_next(1:3)';
    TSX_oe_hist(idx,:)   = TSX_oe_next;
    TDX_oe_hist(idx,:)   = TDX_oe_next;
    rel_oe_hist(idx,:)   = TSX_oe_next(1) * compute_roes(TSX_oe_next, TDX_oe_next)';
    
    % 7) Prepare for next iteration
    TSX_ECI_cur = TSX_ECI_next;
    TDX_ECI_cur = TDX_ECI_next;
end

%% Plotting results
% RTN projections
guidata = rel_oe_hist;  % rename for clarity
figure; subplot(1,3,1); plot(guidata(:,2), guidata(:,1)); axis equal; grid on;
xlabel('r_T [m]'); ylabel('r_R [m]'); title('TR plane');
subplot(1,3,2); plot(guidata(:,3), guidata(:,1)); axis equal; grid on;
xlabel('r_N [m]'); ylabel('r_R [m]'); title('NR plane');
subplot(1,3,3); plot(guidata(:,2), guidata(:,3)); axis equal; grid on;
xlabel('r_T [m]'); ylabel('r_N [m]'); title('TN plane');

% Relative OE vs orbit
figure; subplot(3,2,1); plot(orbit_num, rel_oe_hist(:,1)); ylabel('a\deltaa [m]'); grid on;
subplot(3,2,2); plot(t_orbit, rel_oe_hist(:,2)); ylabel('a\delta\lambda [m]'); grid on;
subplot(3,2,3); plot(t_orbit, rel_oe_hist(:,3)); ylabel('a\deltae_x [m]'); grid on;
subplot(3,2,4); plot(t_orbit, rel_oe_hist(:,4)); ylabel('a\deltae_y [m]'); grid on;
subplot(3,2,5); plot(t_orbit, rel_oe_hist(:,5)); ylabel('a\deltai_x [m]'); grid on;
subplot(3,2,6); plot(t_orbit, rel_oe_hist(:,6)); ylabel('a\deltai_y [m]'); grid on;

% Cumulative Δv budget
cum_dv = cumsum(sum(abs(dv_hist),1));
figure; plot(orbit_num, cum_dv); grid on;
xlabel('Orbit Number'); ylabel('Cumulative Δv [m/s]');
total_dv = cum_dv(end);
fprintf('Total Δv = %.3f m/s\n', total_dv);
