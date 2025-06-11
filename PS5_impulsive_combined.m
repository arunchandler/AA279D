clc; clear; close all;
addpath mean_osc          % your helper folder
format long g

%% ---------- constants & helpers ----------
global tol Re J2 mu s_d
tol = 1e-10;
Re  = 6378137;             J2  = 1.082626e-3;
mu  = 3.986004418e14;      s_d = 86400;

vcol = @(x) x(:);          % force column
vrow = @(x) x(:).';        % force row
wrap = @(x) wrapTo2Pi(mod(x,2*pi));

%% ---------- chief initial orbit ----------
a_TSX    = 6886536.686;    e_TSX = 0.0001264;
i_TSX    = deg2rad(97.4453);
RAAN_TSX = deg2rad(351.0108);
omega_TSX= deg2rad(101.2452);
M_TSX    = deg2rad(11.6520);

TSX_init = [a_TSX e_TSX i_TSX RAAN_TSX omega_TSX M_TSX];

%% ---------- mission scenario ----------
scenario = input( ...
    ['Select mission scenario:\n' ...
     '  1: Formation keeping only\n' ...
     '  2: Single reconfiguration (M-D2 to M-D3)\n' ...
     '  3: Multiple reconfigurations with keeping\n' ...
     'Your choice: '] );

% Define mission phases
switch scenario
    case 1
        % Formation keeping only
        phases = struct('type', 'keeping', ...
                       'start_time', 0, ...
                       'duration', 2*s_d, ...
                       'qns', [0, 0, 0, 300, 0, 400]);
        sim_days = 2;
        
    case 2
        % Single reconfiguration
        phases = [struct('type', 'keeping', ...
                        'start_time', 0, ...
                        'duration', 0.5*s_d, ...
                        'qns', [0, 0, 0, 300, 0, 400]);
                 struct('type', 'reconfig', ...
                        'start_time', 0.5*s_d, ...
                        'duration', 0.1*s_d, ...
                        'qns', [0, 0, 0, 300, 0, 500])];
        sim_days = 1;
        
    case 3
        % Multiple reconfigurations
        phases = [struct('type', 'keeping', ...
                        'start_time', 0, ...
                        'duration', 1*s_d, ...
                        'qns', [0, 0, 0, 300, 0, 400]);
                 struct('type', 'reconfig', ...
                        'start_time', 1*s_d, ...
                        'duration', 0.1*s_d, ...
                        'qns', [0, 0, 0, 300, 0, 500]);
                 struct('type', 'keeping', ...
                        'start_time', 1.1*s_d, ...
                        'duration', 0.9*s_d, ...
                        'qns', [0, 0, 0, 300, 0, 500]);
                 struct('type', 'reconfig', ...
                        'start_time', 2*s_d, ...
                        'duration', 0.1*s_d, ...
                        'qns', [0, 0, 0, 500, 0, 300]);
                 struct('type', 'keeping', ...
                        'start_time', 2.1*s_d, ...
                        'duration', 0.9*s_d, ...
                        'qns', [0, 0, 0, 500, 0, 300])];
        sim_days = 3;
end

%% ---------- initialize deputy ----------
initial_qns = phases(1).qns;
[TDX_init,~] = safe_qns2oe(TSX_init, initial_qns');

%% ---------- formation keeping parameters ----------
de_max_m = 10;             % Maximum eccentricity error [m] - back to original values
di_max_m = 2;              % Maximum inclination error [m] - back to original values
de_max = de_max_m/a_TSX;
di_max = di_max_m/a_TSX;
cooldown = 2*2*pi/sqrt(mu/a_TSX^3);  % 2-orbit cooldown - back to original

%% ---------- grids ----------
n = sqrt(mu/a_TSX^3);
T = 2*pi/n;
t_grid = linspace(0, sim_days*s_d, 1e4+1).';
dt = t_grid(2) - t_grid(1);

%% ---------- state arrays ----------
N = numel(t_grid);
TSX = zeros(N,6);    TSX(1,:) = TSX_init;
TDX = zeros(N,6);    TDX(1,:) = TDX_init;
rel = zeros(N,6);
rel(1,:) = a_TSX * compute_roes(TSX_init, TDX_init)';

%% ---------- bookkeeping ----------
pending = [];              % Burn queue
burn_t = [];              % Executed burn times
burn_dv = [];             % Executed burn magnitudes
burn_vecs = [];           % Executed burn components
burn_types = {};          % Burn types (keeping/reconfig)
cum_dv = 0;
maneuver_log = [];
last_plan = -Inf;
current_phase_idx = 1;
reconfig_done = false;

%% ---------- MAIN LOOP ----------
for k = 2:N
    t_prev = t_grid(k-1);
    t_cur = t_grid(k);
    
    % Propagate chief and deputy
    [~, c_tmp] = ode4(@compute_rates_GVE_J2, [t_prev t_cur]', TSX(k-1,:)', dt);
    [~, d_tmp] = ode4(@compute_rates_GVE_J2, [t_prev t_cur]', TDX(k-1,:)', dt);
    chief = vrow(c_tmp(end,:));
    dep = vrow(d_tmp(end,:));
    
    % Execute queued burns
    while ~isempty(pending) && pending(1).t <= t_cur + 1e-9
        burn = pending(1);
        
        % Rewind to burn time
        [~, c_tmp] = ode4(@compute_rates_GVE_J2, [t_prev burn.t]', chief', dt/10);
        [~, d_tmp] = ode4(@compute_rates_GVE_J2, [t_prev burn.t]', dep', dt/10);
        chief = vrow(c_tmp(end,:));
        dep = vrow(d_tmp(end,:));
        
        % Apply burn
        rv_dep = oe2rv(dep, mu);
        r_dep = rv_dep(1:3);
        v_dep = rv_dep(4:6);
        dv_eci = dv_RTN2ECI(r_dep, v_dep, burn.dv_rtn);
        dep = vrow(kick(dep, dv_eci));
        
        % Record burn
        burn_t = [burn_t; burn.t];
        burn_dv = [burn_dv; norm(dv_eci)];
        burn_vecs = [burn_vecs; get_rtn_components(dep, dv_eci)];
        burn_types{end+1} = burn.type;
        cum_dv = cum_dv + norm(dv_eci);
        
        fprintf('Executed %s burn (%s) at t=%.1f s, |Δv|=%.2f mm/s\n', ...
                burn.tag, burn.type, burn.t, norm(dv_eci)*1e3);
        
        % Propagate to current time
        [~, c_tmp] = ode4(@compute_rates_GVE_J2, [burn.t t_cur]', chief', dt/10);
        [~, d_tmp] = ode4(@compute_rates_GVE_J2, [burn.t t_cur]', dep', dt/10);
        chief = vrow(c_tmp(end,:));
        dep = vrow(d_tmp(end,:));
        
        pending(1) = [];
        t_prev = burn.t;
    end
    
    % Store state
    TSX(k,:) = chief;
    TSX(k,6) = wrapTo2Pi(TSX(k,6));
    TDX(k,:) = dep;
    TDX(k,6) = wrapTo2Pi(TDX(k,6));
    rel(k,:) = a_TSX * compute_roes(chief, dep)';
    
    % Determine current phase
    for p = 1:length(phases)
        if t_cur >= phases(p).start_time && ...
           t_cur < phases(p).start_time + phases(p).duration
            current_phase_idx = p;
            break;
        end
    end
    
    if current_phase_idx > length(phases)
        continue;
    end
    
    current_phase = phases(current_phase_idx);
    
    % Plan maneuvers based on phase type
    if isempty(pending) && (t_cur - last_plan) >= cooldown
        
        % Get nominal values for current phase
        qns_nom = current_phase.qns';
        e_nom_hat = qns_nom(3:4) / a_TSX;
        i_nom_hat = qns_nom(5:6) / a_TSX;
        
        if strcmp(current_phase.type, 'keeping')
            % Plan formation keeping maneuvers
            [new_burns, new_log] = plan_formation_keeping(chief, dep, ...
                e_nom_hat, i_nom_hat, de_max, di_max, t_cur, a_TSX);
            
            if ~isempty(new_burns)
                pending = [pending; new_burns];
                maneuver_log = [maneuver_log; new_log];
                last_plan = t_cur;
                
                fprintf('\n=== Formation Keeping at t=%.1f s ===\n', t_cur);
                fprintf('Planned %d burns\n', length(new_burns));
            end
            
        elseif strcmp(current_phase.type, 'reconfig')
            % Check if reconfiguration for this phase is done
            phase_start = current_phase.start_time;
            already_done = false;
            for i = 1:length(burn_t)
                if burn_t(i) >= phase_start && strcmp(burn_types{i}, 'reconfig')
                    already_done = true;
                    break;
                end
            end
            
            if ~already_done
                % Plan reconfiguration maneuvers
                roe_target_m = qns_nom;
                [new_burns] = plan_reconfiguration(chief, dep, ...
                    roe_target_m, t_cur, a_TSX);
                
                if ~isempty(new_burns)
                    pending = [pending; new_burns];
                    last_plan = t_cur;
                    
                    fprintf('\n=== Reconfiguration at t=%.1f s ===\n', t_cur);
                    fprintf('Target: [%g %g %g %g %g %g] m\n', roe_target_m);
                    fprintf('Planned %d burns\n', length(new_burns));
                end
            end
        end
        
        % Sort pending burns by time
        if ~isempty(pending)
            [~, idx] = sort([pending.t]);
            pending = pending(idx);
        end
    end
end

%% ---------- RESULTS ----------
fprintf('\n=== MISSION SUMMARY ===\n');
fprintf('Total burns: %d\n', length(burn_t));
fprintf('Total Δv = %.3f m/s (%.3f m/s per day)\n', cum_dv, cum_dv/sim_days);

% Count burns by type
n_keeping = sum(strcmp(burn_types, 'keeping'));
n_reconfig = sum(strcmp(burn_types, 'reconfig'));
fprintf('Formation keeping burns: %d\n', n_keeping);
fprintf('Reconfiguration burns: %d\n', n_reconfig);

%% ---------- PLOTS ----------

% ROE evolution
figure('Name', 'ROE Evolution');
labels = {'a·δa [m]','a·δλ [m]','a·δe_x [m]','a·δe_y [m]',...
          'a·δi_x [m]','a·δi_y [m]'};
colors = lines(length(phases));

for i = 1:6
    subplot(3,2,i)
    plot(t_grid/3600, rel(:,i), 'LineWidth', 1.2);
    hold on
    
    % Plot phase boundaries and targets
    for p = 1:length(phases)
        t_start = phases(p).start_time / 3600;
        t_end = (phases(p).start_time + phases(p).duration) / 3600;
        yval = phases(p).qns(i);
        
        plot([t_start t_end], [yval yval], '--', ...
             'Color', colors(p,:), 'LineWidth', 1.5);
        
        if p < length(phases)
            plot([t_end t_end], [yval phases(p+1).qns(i)], 'k:', 'LineWidth', 1);
        end
    end
    
    ylabel(labels{i}); 
    grid on
    if i == 1
        title('Relative Orbital Elements vs Time');
    end
    if i > 4
        xlabel('Time [hours]');
    end
end

% Maneuver timeline
figure('Name', 'Maneuver Timeline');
subplot(3,1,1)
idx_keep = strcmp(burn_types, 'keeping');
idx_reconfig = strcmp(burn_types, 'reconfig');
stem(burn_t(idx_keep)/3600, burn_vecs(idx_keep,1)*1000, 'b', 'LineWidth', 1.5);
hold on
stem(burn_t(idx_reconfig)/3600, burn_vecs(idx_reconfig,1)*1000, 'r', 'LineWidth', 1.5);
ylabel('δv_R [mm/s]'); 
grid on
title('Radial Burns');
legend('Keeping', 'Reconfig', 'Location', 'best');

subplot(3,1,2)
stem(burn_t(idx_keep)/3600, burn_vecs(idx_keep,2)*1000, 'b', 'LineWidth', 1.5);
hold on
stem(burn_t(idx_reconfig)/3600, burn_vecs(idx_reconfig,2)*1000, 'r', 'LineWidth', 1.5);
ylabel('δv_T [mm/s]'); 
grid on
title('Tangential Burns');

subplot(3,1,3)
stem(burn_t(idx_keep)/3600, burn_vecs(idx_keep,3)*1000, 'b', 'LineWidth', 1.5);
hold on
stem(burn_t(idx_reconfig)/3600, burn_vecs(idx_reconfig,3)*1000, 'r', 'LineWidth', 1.5);
ylabel('δv_N [mm/s]'); 
xlabel('Time [hours]');
grid on
title('Normal Burns');

% Cumulative Δv
figure('Name', 'Cumulative Δv');
cum_hist = zeros(N,1);
for i = 1:N
    cum_hist(i) = sum(burn_dv(burn_t <= t_grid(i)));
end
plot(t_grid/3600, cum_hist*1000, 'LineWidth', 1.6);
grid on
ylabel('Cumulative Δv [mm/s]');
xlabel('Time [hours]');
title('Total Δv Consumption');

% Add phase boundaries
hold on
for p = 2:length(phases)
    t_phase = phases(p).start_time / 3600;
    plot([t_phase t_phase], ylim, 'k--', 'LineWidth', 1);
end

function oeNew = kick(oeRow, dvECI)
    global mu
    
    % Helper functions
    vrow = @(x) x(:).';        % force row
    
    rv = oe2rv(oeRow, mu);
    rv(4:6) = rv(4:6) + dvECI(:);
    oe = safe_rv2oe(rv, mu);
    oeNew = [oe(1:5) safe_true2mean(oe(6), oe(2))];
end

function rtn_components = get_rtn_components(oe_dep, dv_vec)
    global mu
    
    % Helper functions
    vrow = @(x) x(:).';        % force row
    
    rv_dep = oe2rv(oe_dep, mu);
    r_dep = rv_dep(1:3);
    v_dep = rv_dep(4:6);
    [Rhat, That, Nhat] = eci2rtn_dir(r_dep, v_dep);
    
    dv_R = dot(dv_vec, Rhat);
    dv_T = dot(dv_vec, That);
    dv_N = dot(dv_vec, Nhat);
    
    rtn_components = [dv_R dv_T dv_N];
end

function [R,T,N] = eci2rtn_dir(r, v)
    % Helper functions
    vrow = @(x) x(:).';        % force row
    
    R = r(:).' / norm(r);
    h = cross(r, v);
    N = h(:).' / norm(h);
    T = cross(N, R);
end

function delta_v_ECI = dv_RTN2ECI(r_ECI, v_ECI, delta_v_RTN)
    % Helper functions
    vcol = @(x) x(:);          % force column
    
    % Normalize R = radial direction (along position)
    R_hat = r_ECI / norm(r_ECI);
    R_hat = R_hat(:);

    % N = orbit normal direction = r × v
    N = cross(r_ECI, v_ECI);
    N_hat = N / norm(N);
    N_hat = N_hat(:);

    T_hat = cross(N_hat, R_hat);
    T_hat = T_hat(:);

    % Construct rotation matrix from RTN to ECI
    Q_RTN2ECI = [R_hat, T_hat, N_hat];
    delta_v_ECI = Q_RTN2ECI * delta_v_RTN;
end

function [burns, maneuver_log] = plan_formation_keeping(chief, dep, e_nom_hat, i_nom_hat, de_max, di_max, t_cur, a_chief)
    % Plan formation keeping maneuvers using tangential and normal burns
    % Following the working PS5_Analytical.m approach
    global Re J2 mu tol
    
    % Helper functions
    vcol = @(x) x(:);          % force column
    vrow = @(x) x(:).';        % force row
    wrap = @(x) wrapTo2Pi(mod(x,2*pi));
    
    burns = [];
    maneuver_log = [];
    
    % Get current ROE (dimensionless)
    roe_hat = compute_roes(chief, dep)';
    de_hat = vcol(roe_hat(3:4));
    di_hat = vcol(roe_hat(5:6));
    
    % Calculate errors
    de_err = norm(de_hat - e_nom_hat);
    di_err = norm(di_hat - i_nom_hat);
    
    % Check if maneuvers are needed
    tol_maneuver = 1e-10;
    need_e_maneuver = de_err > (de_max + tol_maneuver);
    need_i_maneuver = di_err > (di_max + tol_maneuver);
    
    if ~need_e_maneuver && ~need_i_maneuver
        return; % No maneuvers needed
    end
    
    % Calculate J2-perturbed mean motion
    n = sqrt(mu/a_chief^3);
    eta = sqrt(1 - chief(2)^2);
    gamma = 0.5*J2*(Re/a_chief)^2/eta^4;
    P = 3*cos(chief(3))^2 - 1;
    Q = 5*cos(chief(3))^2 - 1;
    kappa = 0.75*J2*Re^2*sqrt(mu)/(a_chief^(7/2)*eta^4);
    nbar = n + kappa*(eta*P + Q);
    
    % Current argument of latitude
    u_c = wrapTo2Pi(mean2true(chief(6), chief(2), tol) + chief(5));
    
    %% Plan eccentricity control (tangential burns) - like PS5_Analytical.m
    if need_e_maneuver
        % Target vectors (following PS5_Analytical.m approach)
        phi_d = 1.5*gamma*(5*cos(chief(3))^2 - 1);
        dphi = sign(phi_d)*asin(de_max/norm(e_nom_hat)) + phi_d*pi;
        Rz = [cos(dphi) -sin(dphi); sin(dphi) cos(dphi)];
        de_man = Rz*e_nom_hat;
        de_req = vcol(de_man) - de_hat;
        
        % Burn epochs
        uM1 = wrapTo2Pi(atan2(de_req(2), de_req(1)));
        uM2 = wrapTo2Pi(uM1 + pi);
        t1 = find_next_burn_time(uM1, u_c, t_cur, nbar);
        t2 = find_next_burn_time(uM2, u_c, t_cur, nbar);
        
        % Δv magnitudes - simplified Damico maneuver (like PS5_Analytical.m)
        dvt1 = n*a_chief/4 * norm(de_req);
        dvt2 = -dvt1;  % Opposite burn for eccentricity control
        
        % Add numerical bounds checking (like PS5_Analytical.m)
        max_dv = 0.1;  % Maximum 10 cm/s per burn
        dvt1 = sign(dvt1) * min(abs(dvt1), max_dv);
        dvt2 = sign(dvt2) * min(abs(dvt2), max_dv);
        
        % Add burns to queue
        burns = [burns;
                struct('t', t1, 'dv_rtn', [0; dvt1; 0], 'tag', 'T1', 'type', 'keeping');
                struct('t', t2, 'dv_rtn', [0; dvt2; 0], 'tag', 'T2', 'type', 'keeping')];
    end
    
    %% Plan inclination control (normal burn) - like PS5_Analytical.m
    if need_i_maneuver
        % Target inclination vector (with J2 drift compensation)
        di_max_new = di_max + abs(3*gamma*di_hat(1)*pi*sin(chief(3))^2);
        di_man = [i_nom_hat(1); i_nom_hat(2) - sign(di_hat(1))*di_max_new];
        di_req = vcol(di_man) - di_hat;
        
        % Burn location
        uN = wrapTo2Pi(atan2(di_req(2), di_req(1)));
        tn = find_next_burn_time(uN, u_c, t_cur, nbar);
        
        % Normal burn magnitude
        dvn = n*a_chief * norm(di_req);
        
        % Add numerical bounds checking to prevent excessive burns
        max_dv_n = 0.5;  % Maximum 50 cm/s per burn for reconfiguration
        dvn = sign(dvn) * min(abs(dvn), max_dv_n);
        
        % Add burn to queue
        burns = [burns;
                struct('t', tn, 'dv_rtn', [0; 0; dvn], 'tag', 'N', 'type', 'keeping')];
    end
    
    % Log the maneuver plan
    dvt1_log = 0; dvt2_log = 0; dvn_log = 0;
    if need_e_maneuver
        dvt1_log = dvt1; dvt2_log = dvt2;
    end
    if need_i_maneuver
        dvn_log = dvn;
    end
    
    maneuver_log = struct('t_det', t_cur, ...
                         'de_error', de_err*a_chief, ...
                         'di_error', di_err*a_chief, ...
                         'dvt1', dvt1_log, ...
                         'dvt2', dvt2_log, ...
                         'dvn', dvn_log);
end

function [burns] = plan_reconfiguration(chief, dep, roe_target_m, t_cur, a_chief)
    % Plan reconfiguration maneuvers using radial and normal burns
    global Re J2 mu tol
    
    % Helper functions
    vcol = @(x) x(:);          % force column
    vrow = @(x) x(:).';        % force row
    wrap = @(x) wrapTo2Pi(mod(x,2*pi));
    
    burns = [];
    
    % Get current ROE in meters
    roe_current_m = a_chief * compute_roes(chief, dep)';
    
    % Convert to dimensionless
    roe_current = vrow(roe_current_m / a_chief);
    roe_target = vrow(roe_target_m / a_chief);
    
    % Required changes
    de_req = vcol(roe_target(3:4) - roe_current(3:4))
    di_req = vcol(roe_target(5:6) - roe_current(5:6))
    dlambda_req = roe_target(2) - roe_current(2);


    assert(isequal(size(de_req),[2 1]), 'de\_req not 2×1');
    assert(isequal(size(di_req),[2 1]), 'di\_req not 2×1');

    fprintf('\n[reconfig-planner]  t = %.0f s\n', t_cur);

    fprintf('   current ROE [m]: %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n', ...
            a_chief*roe_current);
    
    fprintf('   target  ROE [m]: %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n', ...
            roe_target_m);
    
    fprintf('   dROE req  [m]: %7.1f %7.1f %7.1f %7.1f %7.1f %7.1f\n', ...
            a_chief*[0 dlambda_req vrow(de_req) vrow(di_req)]);
            
    % Calculate J2-perturbed mean motion
    n = sqrt(mu/a_chief^3);
    eta = sqrt(1 - chief(2)^2);
    P = 3*cos(chief(3))^2 - 1;
    Q = 5*cos(chief(3))^2 - 1;
    kappa = 0.75*J2*Re^2*sqrt(mu)/(a_chief^(7/2)*eta^4);
    nbar = n + kappa*(eta*P + Q);
    
    % Current argument of latitude
    u_c = wrapTo2Pi(mean2true(chief(6), chief(2), tol) + chief(5));

    if dot([cos(u_c); sin(u_c)], de_req) < 0
        de_req = -de_req;                    % ← FIX: use the opposite branch
    end

    %% Plan in-plane reconfiguration (radial burns)
    if norm(de_req) > 1e-6 || abs(dlambda_req) > 1e-6
        % Burn locations (reference implementation convention)
        xi = atan2(de_req(1), de_req(2));
        
        % Find next available burn times
        k = 0;
        t1 = -100;
        t2 = -100;
        while t1 < t_cur || t2 < t_cur
            uM1 = xi + k*pi;
            uM2 = uM1 + pi;
            t1 = t_cur + mod(uM1 - u_c, 2*pi)/nbar;
            t2 = t_cur + mod(uM2 - u_c, 2*pi)/nbar;
            k = k + 1;
            if k > 10
                error('Could not find valid burn times');
            end
        end
        
        uM1 = wrapTo2Pi(uM1);
        uM2 = wrapTo2Pi(uM2);
        
        % === NEW: compensate J2-induced drift of mean argument of latitude ===
        Delta_t = t2 - t1;
        gamma   = 0.5 * J2 * (Re/a_chief)^2 / eta^4;
        du_J2   = -12 * gamma * sin(2*chief(3)) * roe_current(5) * nbar * Delta_t;
        delta_u_tot = dlambda_req - du_J2;     % net change we *really* need
        dlambda_req = delta_u_tot;             % overwrite the old value
        % Radial burn magnitudes
        de_norm = norm(de_req);
        dvr1 = n*a_chief/2 * (-dlambda_req/2 + de_norm);
        dvr2 = n*a_chief/2 * (-dlambda_req/2 - de_norm);
        
        % Add numerical bounds checking to prevent excessive burns
        % max_dv = 0.5;  % Maximum 50 cm/s per burn for reconfiguration
        % dvr1 = sign(dvr1) * min(abs(dvr1), max_dv);
        % dvr2 = sign(dvr2) * min(abs(dvr2), max_dv);

        % =========== DEBUG block ② ==============================================
        fprintf('   R1 at %.0f s  Δv_R = %+6.2f cm/s\n', t1, dvr1*100);
        fprintf('   R2 at %.0f s  Δv_R = %+6.2f cm/s\n', t2, dvr2*100);
        % ========================================================================
        
        % Add burns to queue
        burns = [burns;
                struct('t', t1, 'dv_rtn', [dvr1; 0; 0], 'tag', 'R1', 'type', 'reconfig');
                struct('t', t2, 'dv_rtn', [dvr2; 0; 0], 'tag', 'R2', 'type', 'reconfig')];
    end
    
    %% -------- OUT-OF-PLANE reconfiguration (single N burn, DEBUG) ---------
    if norm(di_req) > 1e-6
        u_c = wrapTo2Pi( chief(6) + chief(5) ); 
        t_ref = isempty(burns)*t_cur + ~isempty(burns)*max([burns.t]);
        % present true argument-of-latitude
        u_now = wrapTo2Pi(mean2true(chief(6),chief(2),tol) + chief(5));
        
        % latitude the chief will have at t_ref
        u_ref = wrapTo2Pi(u_now + nbar*(t_ref - t_cur));
        % latitude that should give +Δi_y  & 0 Δi_x
        uN_base = wrapTo2Pi( atan2(di_req(2), di_req(1)));
    
        % --- DEBUG ①: geometry check ---------------------------------------
        fprintf('   (DEBUG) u_c      = %+7.4f rad\n', u_c);
        fprintf('   (DEBUG) uN_base  = %+7.4f rad\n', uN_base);
        fprintf('   (DEBUG) dir Gauss= [%6.3f  %6.3f]\n', ...
                -cos(uN_base), -sin(uN_base));
        fprintf('   (DEBUG) Δi_req n = [%6.3f  %6.3f]\n', ...
                di_req/norm(di_req));
        % -------------------------------------------------------------------
    
        k = 0;
        while true
            tn = t_ref + mod(uN_base + 2*pi*k - u_ref , 2*pi)/nbar;
            if tn > t_ref, break; end
            k = k + 1;
        end
    
        dvn = n*a_chief*norm(di_req);
    
        % --- DEBUG ②: epoch & predicted update -----------------------------
        dt_sec = tn - t_ref;
        fprintf('   (DEBUG) tn-t_ref = %.1f s  (%.2f rev)\n', dt_sec, dt_sec*nbar/(2*pi));
        fprintf('   (DEBUG) Δδi_pred = [%7.3f  %7.3f] m\n', ...
                a_chief*[-cos(uN_base); -sin(uN_base)]*dvn/(n*a_chief));
        % -------------------------------------------------------------------
    
        burns = [burns;
                 struct('t', tn, 'dv_rtn', [0;0;dvn], ...
                        'tag','N','type','reconfig')];
    end
end

function t_burn = find_next_burn_time(u_tgt, u_now, t_now, nbar)
    % Find the next time when the satellite reaches the target argument of latitude
    
    % Helper functions
    wrap = @(x) wrapTo2Pi(mod(x,2*pi));
    
    k = 0;
    t_burn = -100;
    
    while t_burn < t_now
        u_target = u_tgt + k*2*pi;
        t_burn = t_now + (u_target - u_now)/nbar;
        k = k + 1;
        
        if k > 100
            error('Could not find valid burn time within 100 revolutions');
        end
    end
end

function [chief_new, dep_new] = execute_burn(chief, dep, burn, t_from, dt_fine)
    global mu

    % propagate to burn epoch -------------------------------------------
    [~, c_tmp] = ode4(@compute_rates_GVE_J2,[t_from burn.t]',chief',dt_fine);
    [~, d_tmp] = ode4(@compute_rates_GVE_J2,[t_from burn.t]', dep', dt_fine);
    chief_at = c_tmp(end,:)';             % row→column
    dep_at   = d_tmp(end,:)';

    % --- DEBUG: where is the chief right now? ---------------------------
    aC  = chief_at(1);  eC = chief_at(2);
    wC  = chief_at(5);
    vC  = mean2true(chief_at(6), eC, 1e-12);   % true anomaly
    uC  = mod(vC + wC , 2*pi);                 % true argument of latitude
    fprintf('      (debug exec)  burn.lat true = %.4f rad  (%.1f°)\n', ...
            uC, uC*180/pi);
    % -------------------------------------------------------------------

    % apply dv -----------------------------------------------------------
    rv  = oe2rv(dep_at',mu);
    dvE = dv_RTN2ECI(rv(1:3),rv(4:6), burn.dv_rtn);
    dep_after = kick(dep_at', dvE)';

    chief_new = chief_at;     dep_new = dep_after;
end
