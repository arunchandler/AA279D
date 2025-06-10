clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day 

use_normal_burn = true;   % << set true to re-enable

vcol = @(x) x(:);      % force column
vrow = @(x) x(:).';    % force row

%% Problem 1 - Real-time Formation Control

% Define target ROE values based on PDF Table 2.2
% aδenom = 500 m, φnom = 80°, aδinom = 300 m, θnom = 50°
roe_target = [0; 0; 500*cos(deg2rad(80)); 500*sin(deg2rad(80)); 300*cos(deg2rad(50)); 300*sin(deg2rad(50))];

scenario = input( ...
    ['Select rel_qns_init scenario:\n' ...
     '  1: DEM formation (e_term=50/√2)\n' ...
     '  2: Example "circular" formation\n' ...
     '  3: Custom (enter manually)\n' ...
     'Your choice: '] );

switch scenario
    case 1
        % DEM - use smaller relative eccentricity
        e_term = 50/sqrt(2);  % Reduced from 122/sqrt(2)
        rel_qns_init = [0, 340, e_term, e_term, 0, 256];
        scenario_name = 'DEM';

    case 2
        % Pursuit - use smaller relative eccentricity
        e_term = 100/sqrt(2);  % Reduced from 4982/sqrt(2)
        rel_qns_init = [0, -76050, e_term, e_term, 80, 0 ];  
        scenario_name = 'Pursuit';

    case 3
        % Large cross-track - use smaller relative eccentricity
        e_term = 200/sqrt(2);  % Reduced from 3600/sqrt(2)
        rel_qns_init = [0, 0, e_term, e_term, 250, 0 ];  
        scenario_name = 'Large cross-track';

    case 4
        % Short baseline - use smaller relative eccentricity
        e_term = 100/sqrt(2);  % Reduced from 250/sqrt(2)
        rel_qns_init = [0, -340, e_term, e_term, 250, 0 ];  
        scenario_name = 'Short baseline';

    case 5
        % Damico separation - use target values as initial
        % Based on PDF Table 2.2: aδenom = 500 m, φnom = 80°, aδinom = 300 m, θnom = 50°
        % Use the same values as roe_target to ensure consistency
        rel_qns_init = roe_target';
        scenario_name = 'Damico Paper Initialization';

    otherwise
        error('Invalid scenario');
end

fprintf('Running scenario "%s": rel_qns_init = [%g %g %g %g %g %g]\n\n', ...
        scenario_name, rel_qns_init);

% Initial orbital elements
a_TSX_init    = 6886536.686;  % m
e_TSX_init    = 0.0001264;
i_TSX_init    = deg2rad(97.4453);
RAAN_TSX_init = deg2rad(351.0108);
omega_TSX_init= deg2rad(101.2452);
M_TSX_init    = deg2rad(11.6520);
nu_TSX_init   = mean2true(M_TSX_init, e_TSX_init, tol);
u_TSX_init    = nu_TSX_init + omega_TSX_init;

TSX_init_oe = [a_TSX_init, e_TSX_init, i_TSX_init, ...
               RAAN_TSX_init, omega_TSX_init, M_TSX_init];

% Map to TDX via chosen rel_qns_init
[TDX_init_oe, qns_success] = safe_qns2oe(TSX_init_oe, rel_qns_init);

% Debug the initial ROE conversion
debug_roe_conversion(TSX_init_oe, TDX_init_oe, rel_qns_init);

if ~qns_success
    fprintf('WARNING: Initial ROE conversion required eccentricity clamping!\n');
    fprintf('Consider reducing the relative eccentricity components in rel_qns_init.\n\n');
end

% Timing parameters
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init^3);
T          = 2*pi/n;
n_orbit = ceil((30 * 86400) / T);  % Number of orbits in 30 days
sim_days = 3;                % §2.4 uses three days
tend     = sim_days*86400;   
num_points = 10000;  % Increase for better resolution
dt = (tend - tstart)/(num_points-1);
t_grid = linspace(tstart, tend, num_points).';
t_orbit = t_grid / T;

% Logging for performance plots
err_hist   = zeros(num_points, 6);     % a⋅(δ – δ_nom)  [m]
e_norm     = zeros(num_points, 1);     % ‖a·δe‖
i_norm     = zeros(num_points, 1);     % ‖a·δi‖

burn_t     = [];                       % time [s] of every burn
burn_dv    = [];                       % |Δv|  per burn  [m/s]
burn_vecs  = [];                       % 3-component Δv  (for direction stats)
cum_dv     = 0;                        % running total
cum_dv_hist = zeros(num_points,1);

%% === REAL-TIME FORMATION CONTROL ===

a_chief = a_TSX_init;                     % chief SMA  [m]
n       = sqrt(mu/a_chief^3);             % mean motion [rad s-1]
T       = 2*pi/n;                         % orbital period [s]
T_half  = T/2;                            % half orbital period [s]

% Control window - use larger values to avoid immediate violations
delta_e_max_m = 2.0;                    % half-window ‖a·δe‖ [m]
delta_i_max_m = 2.0;                    % half-window ‖a·δi‖ [m]

% Convert to dimensionless
e_nom_hat = roe_target(3:4)/a_chief;      % δe_nom   (dimensionless)
i_nom_hat = roe_target(5:6)/a_chief;      % δi_nom   (dimensionless)
delta_u_nom = - i_nom_hat(2) / tan(i_TSX_init);
de_max    = delta_e_max_m/a_chief;        % ‖δe‖ window (dimensionless)
di_max    = delta_i_max_m/a_chief;        % ‖δi‖ window (dimensionless)

% DEBUG: Print initial ROE values to understand the drift
fprintf('\n=== INITIAL ROE ANALYSIS ===\n');
roe_initial = compute_roes(TSX_init_oe, TDX_init_oe)';
roe_initial_m = a_chief * roe_initial;  % Convert to meters
fprintf('Initial ROE (m): [%.2f %.2f %.2f %.2f %.2f %.2f]\n', roe_initial_m);
fprintf('Target ROE (m):  [%.2f %.2f %.2f %.2f %.2f %.2f]\n', roe_target');

% Initialize real-time control variables
maneuver_log = [];                        % save maneuver data
maneuver_triggered = false;
last_maneuver_time = 0;                   % track last maneuver time
maneuver_cooldown = T;                    % minimum time between maneuvers

% Initialize orbital element arrays
TSX_oe = zeros(num_points, 6);
TDX_oe = zeros(num_points, 6);
rel_oe = zeros(num_points, 6);
rtn    = zeros(num_points, 6);

% Set initial conditions
TSX_oe(1,:) = TSX_init_oe;
TDX_oe(1,:) = TDX_init_oe;

% Real-time control loop with propagation
idx = 1;

while idx <= num_points
    
    % Propagate both satellites with J2
    if idx == 1
        % Initial conditions already set
    else
        % Propagate from previous state
        [~, TSX_temp] = ode4(@compute_rates_GVE_J2, [t_grid(idx-1), t_grid(idx)]', TSX_oe(idx-1,:)', dt);
        [~, TDX_temp] = ode4(@compute_rates_GVE_J2, [t_grid(idx-1), t_grid(idx)]', TDX_oe(idx-1,:)', dt);
        
        TSX_oe(idx,:) = TSX_temp(end,:);
        TDX_oe(idx,:) = TDX_temp(end,:);
        
        % Wrap mean anomaly
        TSX_oe(idx,6) = wrapTo2Pi(TSX_oe(idx,6));
        TDX_oe(idx,6) = wrapTo2Pi(TDX_oe(idx,6));
    end
    
    % Compute relative OEs and RTN
    a1 = TSX_oe(idx,1); e1 = TSX_oe(idx,2); i1 = TSX_oe(idx,3);
    RAAN1 = TSX_oe(idx,4); omega1 = TSX_oe(idx,5);
    M1 = wrapTo2Pi(TSX_oe(idx,6)); u1 = M1 + omega1;

    a2 = TDX_oe(idx,1); e2 = TDX_oe(idx,2); i2 = TDX_oe(idx,3);
    RAAN2 = TDX_oe(idx,4); omega2 = TDX_oe(idx,5);
    M2 = wrapTo2Pi(TDX_oe(idx,6)); u2 = M2 + omega2;

    rel_oe(idx,:) = a1 * compute_roes([a1, e1, i1, RAAN1, omega1, M1], ...
                                      [a2, e2, i2, RAAN2, omega2, M2])';

    r1 = oe2rv(TSX_oe(idx,:), mu);
    r2 = oe2rv(TDX_oe(idx,:), mu);
    rtn(idx,:) = eci2rtn(r1, r2)';
    
    % Current relative elements - IMPORTANT: compute_roes returns DIMENSIONLESS ROEs
    roe_hat  = compute_roes(TSX_oe(idx,:), TDX_oe(idx,:))';  % Already dimensionless!
    da_hat   = roe_hat(1);                % δa (dimensionless)
    dl_hat   = roe_hat(2);                % δλ (dimensionless)
    de_hat   = vcol(roe_hat(3:4));        % δe  (x,y) - ensure column
    di_hat   = vcol(roe_hat(5:6));        % δi  (x,y) - ensure column
    
    % Chief argument of latitude (for timing)
    M_c  = wrapTo2Pi(TSX_oe(idx,6));
    nu_c = mean2true(M_c, TSX_oe(idx,2), tol);
    u_c  = wrapTo2Pi(nu_c + TSX_oe(idx,5));

    % ── rotate nominal vector every sample ─────────────────────────────
    eta = sqrt(1 - e1^2);
    gamma   = 0.5*J2*(Re/a_chief)^2*(1/(eta^4));
    phi_dot = 1.5*gamma*(5*cos(TSX_oe(idx,3))^2 - 1);
    % For i-vector: The drift is in RAAN, not in the i-vector directly!
    % From Eq. 2.30: δΩ̇ = -3γ cos(i) δi_x
    Omega_dot = -3*gamma*cos(TSX_oe(idx,3)) * di_hat(1);  % RAAN drift rate
    theta_dot = -1.5*gamma*sin(2*TSX_oe(idx,3));  % i-vector drift rate
    theta_nom = theta_dot * t_grid(idx);

    % updated phi calculation
    dphi = sign(phi_dot) * asin(de_max/norm(e_nom_hat));
    
    % Total rotation/drift since t=0
    phi_nom = phi_dot * t_grid(idx);
    Omega_drift = Omega_dot * t_grid(idx);
    
    % Rotate e-vector
    Rz_phi = [cos(dphi) -sin(dphi);
              sin(dphi)  cos(dphi)];
    delta_e_man = Rz_phi * vcol(e_nom_hat);

    % i-vector doesn't rotate in (i_x, i_y) space, it stays fixed!
    % The drift appears in the relative RAAN, not in the i-vector components
    i_nom_rot_hat = vcol(i_nom_hat);  % NO DRIFT in i_x, i_y components!

    di_max_new = di_max + abs(3*gamma*di_hat(1)*pi*sin(i1)^2);

    delta_i_man_y = i_nom_hat(2) - sign(di_hat(1))*di_max_new;

    delta_i_man = [i_nom_hat(1); delta_i_man_y];
    
    % Create rotated nominal ROE vector for logging
    %roe_nom_rot_m = [0; 0; a_chief*e_nom_rot_hat; a_chief*i_nom_rot_hat];
    
    % Update error history (before any maneuvers)
    roe_current_m = a_chief * roe_hat;
    err_hist(idx,:) = vrow(roe_current_m) - vrow(roe_target);
    e_norm(idx) = norm(roe_current_m(3:4) - roe_target(3:4));
    i_norm(idx) = norm(roe_current_m(5:6) - roe_target(5:6));

    % 7. ADD DEBUG OUTPUT to track what's happening
    if mod(idx, 1000) == 0  % Every 1000 points
        fprintf('Debug at orbit %.1f:\n', t_grid(idx)/T);
        fprintf('  e_man: [%.2f, %.2f] m\n', a_chief*delta_e_man(1), a_chief*delta_e_man(2));
        fprintf('  i_man: [%.2f, %.2f] m (static)\n', a_chief*delta_i_man(1), a_chief*delta_i_man(2));
        fprintf('  Current δe: [%.2f, %.2f] m\n', a_chief*de_hat(1), a_chief*de_hat(2));
        fprintf('  Current δi: [%.2f, %.2f] m\n', a_chief*di_hat(1), a_chief*di_hat(2));
        fprintf('  phi rotation: %.2f deg\n', rad2deg(dphi));
    end
    
    % WINDOW CHECK using rotated nominal values
    de_error = norm(de_hat - e_nom_hat);  % dimensionless
    di_error = norm(di_hat - i_nom_hat);  % dimensionless
    
    violation = (de_error > de_max) || (di_error > di_max);
    time_since_last_maneuver = t_grid(idx) - last_maneuver_time;

    % Check if we need to trigger a maneuver (with cooldown)
    if violation && ~maneuver_triggered && time_since_last_maneuver >= maneuver_cooldown
        
        fprintf('\n[%6.2f orbits] Control window violation at idx=%d\n', t_grid(idx)/T, idx);
        fprintf('  δe error: %.2f m (limit: %.2f m)\n', de_error*a_chief, delta_e_max_m);
        fprintf('  δi error: %.2f m (limit: %.2f m)\n', di_error*a_chief, delta_i_max_m);
        
        % ===============================================================
        % Post-maneuver targets (Eqs. 2.58-2.60)
        % ===============================================================
        
        % For e-vector: rotate to opposite side of window
        % if norm(e_nom_rot_hat) > 0
        %     dphi = sign(phi_dot) * asin(de_max/norm(e_nom_rot_hat)) + phi_dot*pi;
        % else
        %     dphi = 0;
        % end
        
        % Rz_dphi = [cos(dphi) -sin(dphi);
        %            sin(dphi)  cos(dphi)];
        % de_man_hat = Rz_dphi * e_nom_rot_hat;  % Target on opposite side of window
        % 
        % % For i-vector: shift to opposite side
        % di_man_hat = vcol(i_nom_hat);
        % if norm(i_nom_rot_hat) > 0 && di_error > di_max
        %     % Move to opposite side of window
        %     scale = (norm(i_nom_hat) - di_max) / norm(i_nom_hat);
        %     di_man_hat = scale * i_nom_hat;
        % end

        % Corrections required
        de_req = vcol(delta_e_man) - vcol(de_hat);  % dimensionless 2x1
        di_req = vcol(delta_i_man) - vcol(di_hat);  % dimensionless 2x1
        
        % ===============================================================
        % Optimal maneuver locations
        % ===============================================================
        
        if norm(de_req) > 0
            xi = atan2(de_req(2), de_req(1));  % Direction of REQUIRED CHANGE
        else
            xi = 0;  % Default if no change needed
        end

        
        % For tangential burns, we burn perpendicular to the change direction
        % This is because tangential burns change eccentricity in the perpendicular direction
        u1 = wrapTo2Pi(xi);   % First tangential burn location
        u2 = wrapTo2Pi(u1 + pi);     % Second tangential burn location (180° apart)
        
        % Calculate times for tangential burns
        du1 = wrapTo2Pi(u1 - u_c);
        du2 = wrapTo2Pi(u2 - u_c);
        t1  = t_grid(idx) + du1/n;
        t2  = t_grid(idx) + du2/n;
        
        % Cross-track burn
        dvn = n*a_chief * norm(di_req);  % Eq. 2.64
        if use_normal_burn && di_error > di_max
            % Only correct the excess error
            di_excess = di_error - di_max;
            dvn = n*a_chief * norm(di_req);  % Only correct the excess!
            % Use:
            if norm(di_req) > 0
                un = atan2(di_req(2), di_req(1));
            else
                un = 0;
            end
            
            
            du_n = wrapTo2Pi(un - u_c);
            tn = t_grid(idx) + du_n/n;
        else
            dvn = 0;
            un = 0;
            tn = NaN;
        end
        
        % ===============================================================
        % Compute δa^man for along-track control (Eq. 2.71)
        % ===============================================================
        Delta_t = T;  % maintenance cycle (1 full orbit)
        
        % Current δu (relative mean argument of latitude)
        %delta_u = u1 - u2;  % dimensionless
        %delta_u_nom = -i_nom_hat(2)/ tan(TSX_oe(idx,3));  % nominal target (dimensionless)
        
        % Secular contributions (Eq. 2.69)
        du_J2 = -12*gamma*sin(2*TSX_oe(idx,3))*di_hat(1)*n*Delta_t;
        du_D  = 0;  % neglecting differential drag
        
        % --- along-track quantities -------------------------------------------
        delta_u     = dl_hat;          %  the current along-track ROE  [rad]
        delta_u_nom = 0;               %  0  for centred formation
        ee_req      = norm(de_req);    %  |δe_req|  dimensionless
        % ----------------------------------------------------------------------
        
        bracket = 3*de_max + da_hat ...
                  - (4/(3*pi))*( delta_u + du_J2);   % Eq 2·71
        den     = 2*n*Delta_t - pi;                          % stays the same
        da_man  = -pi/den * bracket;
        
        % ===============================================================
        % Δv magnitudes (Eqs. 2.63-2.64)
        % ===============================================================
        dvt1 = (n*a_chief/4) * ((da_man - da_hat) + norm(de_req));
        dvt2 = (n*a_chief/4) * ((da_man - da_hat) - norm(de_req));

        if de_error > de_max

            dvt1 = (n*a_chief/4) * ((da_man - da_hat) + norm(de_req));
            dvt2 = (n*a_chief/4) * ((da_man - da_hat) - norm(de_req));

            if norm(de_req) > 0
                xi = atan2(de_req(2), de_req(1));
                u1 = wrapTo2Pi(xi);   % First tangential burn location
                u2 = wrapTo2Pi(u1 + pi);     % Second tangential burn location (180° apart)
            else
                u1 = 0;
                u2 = 0;
            end
            
            
            % Calculate times for tangential burns
            du1 = wrapTo2Pi(u1 - u_c);
            du2 = wrapTo2Pi(u2 - u_c);
            t1  = t_grid(idx) + du1/n;
            t2  = t_grid(idx) + du2/n;
        else
            dvt1 = 0;
            dvt2 = 0;
            u1 = 0;
            u2 = 0;
            t1 = NaN;
            t2 = NaN;
        end
        
        % Add safeguards
        max_dv = 10.0;  % Maximum allowed Δv in m/s
        if abs(dvt1) > max_dv
            fprintf('WARNING: Δv_T1 = %.4f m/s exceeds limit. Clamping.\n', dvt1);
            dvt1 = sign(dvt1) * max_dv;
        end
        if abs(dvt2) > max_dv
            fprintf('WARNING: Δv_T2 = %.4f m/s exceeds limit. Clamping.\n', dvt2);
            dvt2 = sign(dvt2) * max_dv;
        end
        
        % Save maneuver data
        maneuver_log = [maneuver_log; struct( ...
            'idx',   idx,        't_det', t_grid(idx), ...
            'dvt1',  dvt1,       'dvt2',  dvt2,        ...
            'dvn',   dvn,        'u1',    u1,          ...
            'u2',    u2,         'un',    un,          ...
            'u_c',   u_c,        'de_error', de_error*a_chief, ...
            'di_error', di_error*a_chief, ...
            'de_man', delta_e_man, 'di_man', delta_i_man, ...
            'da_man', da_man )];
        
        fprintf('  Δv_T1 = %.4f m/s @ u1 = %.2f deg\n', dvt1, rad2deg(u1));
        fprintf('  Δv_T2 = %.4f m/s @ u2 = %.2f deg\n', dvt2, rad2deg(u2));
        if use_normal_burn && dvn > 1e-6
            fprintf('  Δv_N  = %.4f m/s @ un = %.2f deg\n', dvn, rad2deg(un));
        end
        fprintf('  Total Δv = %.4f m/s\n\n', abs(dvt1) + abs(dvt2) + dvn);

        % Replace your maneuver application section with this debugged version:

        % ========= APPLY BURNS ===========================================
        
        % Current states at detection time
        chief_now = TSX_oe(idx,:);
        dep_now   = TDX_oe(idx,:);
        
        % Store the state at detection time for propagation
        chief_det = chief_now;
        dep_det = dep_now;
        t_det = t_grid(idx);
        
        % DEBUG: Print state before maneuvers
        fprintf('DEBUG: Pre-maneuver state at idx=%d:\n', idx);
        roe_before = compute_roes(chief_now, dep_now)';
        fprintf('  ROE before: [%.2f %.2f %.2f %.2f %.2f %.2f] m\n', a_chief*roe_before);
        
        % Determine burn order and apply them
        burn_times = [];
        burn_types = {};  % Track burn types for debugging
        
        if use_normal_burn && dvn > 1e-6 && ~isnan(tn)
            burn_times = [burn_times; tn];
            burn_types{end+1} = 'N';
        end
        burn_times = [burn_times; t1; t2];
        burn_types{end+1} = 'T1';
        burn_types{end+2} = 'T2';
        
        % Sort burns by time
        [burn_times, sort_idx] = sort(burn_times);
        burn_types = burn_types(sort_idx);
        
        % Apply burns in chronological order
        fprintf('DEBUG: Applying %d burns\n', length(burn_times));
        
        for burn_idx = 1:length(burn_times)
            t_burn = burn_times(burn_idx);
            burn_type = burn_types{burn_idx};
            
            % Propagate from current position to burn time
            if burn_idx == 1
                t_from = t_det;
            else
                t_from = burn_times(burn_idx-1);
            end
            
            fprintf('  Burn %d (%s): propagating from t=%.1f to t=%.1f\n', ...
                    burn_idx, burn_type, t_from, t_burn);
            
            [~, chief_prop] = ode4(@compute_rates_GVE_J2, [t_from, t_burn]', chief_now', dt/10);
            [~, dep_prop] = ode4(@compute_rates_GVE_J2, [t_from, t_burn]', dep_now', dt/10);
            chief_burn = chief_prop(end,:);
            dep_burn = dep_prop(end,:);
            
            if strcmp(burn_type, 'N')
                % Normal burn
                dvn_vec = dvn * rtn_axis(chief_burn, dep_burn, 'N');  % Now returns column vector
                dep_burn_new = kick(dep_burn, dvn_vec);
                
                fprintf('    Applied normal burn: Δv_N = %.4f m/s\n', dvn);
                fprintf('    Deputy e before: %.6f, after: %.6f\n', dep_burn(2), dep_burn_new(2));
                
                dep_burn = dep_burn_new;
                
                burn_t = [burn_t; t_burn];
                burn_dv = [burn_dv; dvn];
                burn_vecs = [burn_vecs; dvn_vec'];  % Transpose to store as row
                cum_dv = cum_dv + dvn;
                
            elseif strcmp(burn_type, 'T1')
                % First tangential burn
                dv1_vec = dvt1 * rtn_axis(chief_burn, dep_burn, 'T');  % Now returns column vector
                dep_burn_new = kick(dep_burn, dv1_vec);
                
                fprintf('    Applied T1 burn: Δv_T1 = %.4f m/s\n', dvt1);
                fprintf('    Deputy a before: %.2f, after: %.2f\n', dep_burn(1), dep_burn_new(1));
                
                dep_burn = dep_burn_new;
                
                burn_t = [burn_t; t_burn];
                burn_dv = [burn_dv; abs(dvt1)];
                burn_vecs = [burn_vecs; dv1_vec'];  % Transpose to store as row
                cum_dv = cum_dv + abs(dvt1);
                
            elseif strcmp(burn_type, 'T2')
                % Second tangential burn
                dv2_vec = dvt2 * rtn_axis(chief_burn, dep_burn, 'T');  % Now returns column vector
                dep_burn_new = kick(dep_burn, dv2_vec);
                
                fprintf('    Applied T2 burn: Δv_T2 = %.4f m/s\n', dvt2);
                fprintf('    Deputy a before: %.2f, after: %.2f\n', dep_burn(1), dep_burn_new(1));
                
                dep_burn = dep_burn_new;
                
                burn_t = [burn_t; t_burn];
                burn_dv = [burn_dv; abs(dvt2)];
                burn_vecs = [burn_vecs; dv2_vec'];  % Transpose to store as row
                cum_dv = cum_dv + abs(dvt2);
            end
            
            % Update current states for next propagation
            chief_now = chief_burn;
            dep_now = dep_burn;
        end
        
        % DEBUG: Print state after all maneuvers
        fprintf('DEBUG: Post-maneuver state:\n');
        roe_after = compute_roes(chief_now, dep_now)';
        fprintf('  ROE after: [%.2f %.2f %.2f %.2f %.2f %.2f] m\n', a_chief*roe_after);
        fprintf('  Change in ROE: [%.2f %.2f %.2f %.2f %.2f %.2f] m\n', ...
                a_chief*(roe_after - roe_before));
        
        % Find the grid index right after the last burn
        k_end = find(t_grid > burn_times(end), 1, 'first');
        if isempty(k_end)
            k_end = num_points;
        end
        
        % CRITICAL: Store the final post-burn states for future propagation
        chief_final = chief_now;
        dep_final = dep_now;
        
        % Now fill in all the grid points from detection to end of burns
        fprintf('DEBUG: Filling grid from idx=%d to k_end=%d\n', idx, k_end);
        
        % Use a flag to track if we should skip ahead
        skip_to_idx = k_end;
        
        % Fill the grid - this is a separate operation that doesn't affect the main loop
        for jj = idx+1:k_end
            t_jj = t_grid(jj);
            
            % Find which burn segment this time belongs to
            burn_before = find(burn_times <= t_jj, 1, 'last');
            
            if isempty(burn_before)
                % Before any burns
                [~, c_tmp] = ode4(@compute_rates_GVE_J2, [t_det, t_jj]', chief_det', dt/10);
                [~, d_tmp] = ode4(@compute_rates_GVE_J2, [t_det, t_jj]', dep_det', dt/10);
            else
                % Need to reconstruct state after burns_before burns
                chief_segment = chief_det;
                dep_segment = dep_det;
                
                % Apply all burns up to burn_before
                for bi = 1:burn_before
                    tb = burn_times(bi);
                    
                    % Propagate to burn
                    if bi == 1
                        t_start = t_det;
                    else
                        t_start = burn_times(bi-1);
                    end
                    
                    [~, c_prop] = ode4(@compute_rates_GVE_J2, [t_start, tb]', chief_segment', dt/10);
                    [~, d_prop] = ode4(@compute_rates_GVE_J2, [t_start, tb]', dep_segment', dt/10);
                    chief_segment = c_prop(end,:);
                    dep_segment = d_prop(end,:);
                    
                    % Apply burn to deputy
                    burn_type = burn_types{bi};
                    if strcmp(burn_type, 'N')
                        dv_vec = dvn * rtn_axis(chief_segment, dep_segment, 'N');
                        dep_segment = kick(dep_segment, dv_vec);
                    elseif strcmp(burn_type, 'T1')
                        dv_vec = dvt1 * rtn_axis(chief_segment, dep_segment, 'T');
                        dep_segment = kick(dep_segment, dv_vec);
                    elseif strcmp(burn_type, 'T2')
                        dv_vec = dvt2 * rtn_axis(chief_segment, dep_segment, 'T');
                        dep_segment = kick(dep_segment, dv_vec);
                    end
                end
                
                % Now propagate from last burn to current time
                t_start = burn_times(burn_before);
                [~, c_tmp] = ode4(@compute_rates_GVE_J2, [t_start, t_jj]', chief_segment', dt/10);
                [~, d_tmp] = ode4(@compute_rates_GVE_J2, [t_start, t_jj]', dep_segment', dt/10);
            end
            
            % Store states
            TSX_oe(jj,:) = c_tmp(end,:);
            TDX_oe(jj,:) = d_tmp(end,:);
            TSX_oe(jj,6) = wrapTo2Pi(TSX_oe(jj,6));
            TDX_oe(jj,6) = wrapTo2Pi(TDX_oe(jj,6));
            
            % Update tracking arrays
            roe_jj = compute_roes(TSX_oe(jj,:), TDX_oe(jj,:))';
            rel_oe(jj,:) = a_chief * roe_jj;
            
            % Calculate error relative to target
            err_hist(jj,:) = a_chief*vrow(roe_jj) - vrow(roe_target/a_chief);
            e_norm(jj) = norm(a_chief*roe_jj(3:4) - roe_target(3:4));
            i_norm(jj) = norm(a_chief*roe_jj(5:6) - roe_target(5:6));
            
            % RTN
            r1 = oe2rv(TSX_oe(jj,:), mu);
            r2 = oe2rv(TDX_oe(jj,:), mu);
            rtn(jj,:) = eci2rtn(r1, r2)';
        end
        
        % CRITICAL: Make sure the final states are stored at k_end
        if k_end <= num_points
            TSX_oe(k_end,:) = chief_final;
            TDX_oe(k_end,:) = dep_final;
            
            % Verify the change persisted
            roe_verify = compute_roes(TSX_oe(k_end,:), TDX_oe(k_end,:))';
            fprintf('DEBUG: ROE at k_end=%d: [%.2f %.2f %.2f %.2f %.2f %.2f] m\n', ...
                    k_end, a_chief*roe_verify);
        end
        
        % Set maneuver triggered flag and update time
        last_maneuver_time = t_det;
        maneuver_triggered = true;
        
        % Set a flag instead of modifying idx directly
        % The main loop will check this flag
        idx = k_end - 1;
    end
    
    % Reset maneuver trigger after cooldown
    if maneuver_triggered && time_since_last_maneuver >= maneuver_cooldown
        maneuver_triggered = false;
    end

    idx = idx + 1;
end

% Unpack RTN for plotting
rR = rtn(:,1);
rT = rtn(:,2);
rN = rtn(:,3);

%% === REAL-TIME CONTROL COMPLETE ===

% Print final statistics
fprintf('\n=== REAL-TIME CONTROL SUMMARY ===\n');
fprintf('Total maneuvers executed: %d\n', length(maneuver_log));
fprintf('Total Δv consumed: %.4f m/s\n', cum_dv);
fprintf('Average Δv per maneuver: %.4f m/s\n', cum_dv/max(1, length(maneuver_log)));

fprintf('\n=== ERROR EVOLUTION ANALYSIS ===\n');
% Find key time points
orbit_1 = find(t_grid >= T, 1);
orbit_5 = find(t_grid >= 5*T, 1);
orbit_10 = find(t_grid >= 10*T, 1);
orbit_15 = find(t_grid >= 15*T, 1);
orbit_20 = find(t_grid >= 20*T, 1);
orbit_25 = find(t_grid >= 25*T, 1);
orbit_30 = find(t_grid >= 30*T, 1);

key_indices = [1, orbit_1, orbit_5, orbit_10, orbit_15, orbit_20, orbit_25, orbit_30];
key_indices = key_indices(key_indices <= num_points);

fprintf('%-12s %-12s %-12s %-12s %-12s\n', 'Orbit', 'de_error', 'di_error', '||de||', '||di||');
fprintf('%-12s %-12s %-12s %-12s %-12s\n', '-----', '--------', '--------', '------', '------');

% Need to recalculate gamma and phi_dot for error calculations
eta = sqrt(1-TSX_oe(1,2)^2);
gamma = 0.5*J2*(Re/a_chief)^2*(1/(eta^4));
phi_dot = 1.5*gamma*(5*cos(TSX_oe(1,3))^2 - 1);  % Use initial inclination

for i = 1:length(key_indices)
    idx = key_indices(i);
    if idx <= num_points && idx > 0
        orbit_num = t_grid(idx) / T;
        
        % Get current ROE
        roe_tmp = compute_roes(TSX_oe(idx,:), TDX_oe(idx,:))';
        
        % Calculate rotated nominal at this time
        phi_idx = phi_dot * t_grid(idx);
        Rz_idx = [cos(phi_idx) -sin(phi_idx); sin(phi_idx) cos(phi_idx)];
        e_rot_idx = Rz_idx * vcol(e_nom_hat);
        
        % i-vector doesn't drift in ROE representation
        i_rot_idx = vcol(i_nom_hat);
        
        % Calculate errors in meters
        de_error_actual = norm(a_chief*roe_tmp(3:4) - a_chief*e_rot_idx);
        di_error_actual = norm(a_chief*roe_tmp(5:6) - a_chief*i_rot_idx);
        
        fprintf('%-12.1f %-12.2f %-12.2f %-12.2f %-12.2f\n', ...
                orbit_num, de_error_actual, di_error_actual, ...
                norm(a_chief*roe_tmp(3:4)), norm(a_chief*roe_tmp(5:6)));
    end
end

% Show maneuver details
if ~isempty(maneuver_log)
    fprintf('\n=== MANEUVER DETAILS ===\n');
    fprintf('%-8s %-12s %-12s %-12s %-12s\n', 'Maneuver', 'Orbit', 'de_error', 'di_error', 'Total_dV');
    fprintf('%-8s %-12s %-12s %-12s %-12s\n', '-------', '-----', '--------', '--------', '--------');
    
    for i = 1:length(maneuver_log)
        det = maneuver_log(i);
        orbit_num = det.t_det / T;
        total_dv = abs(det.dvt1) + abs(det.dvt2);
        if isfield(det, 'dvn')
            total_dv = total_dv + det.dvn;
        end
        fprintf('%-8d %-12.2f %-12.2f %-12.2f %-12.4f\n', ...
                i, orbit_num, det.de_error, det.di_error, total_dv);
    end
end

%% ───────────────────────────────────────────────────────────
%%  RESULTS & VISUALISATION (FIXED)
%% ───────────────────────────────────────────────────────────

% Fix the error history - it should be computed for ALL points, not just during maneuvers
% Recompute error history for smooth plotting
for idx = 1:num_points
    % Get current ROE
    roe_current = compute_roes(TSX_oe(idx,:), TDX_oe(idx,:))';
    roe_current_m = a_chief * roe_current;
    
    % Calculate rotated nominal at this time
    eta = sqrt(1 - TSX_oe(idx,2)^2);
    gamma = 0.5*J2*(Re/a_chief)^2*(1/(eta^4));
    phi_dot = 1.5*gamma*(5*cos(TSX_oe(idx,3))^2 - 1);
    phi_nom = phi_dot * t_grid(idx);
    
    % Rotate e-vector nominal
    %dphi = 
    Rz_phi = [cos(phi_nom) -sin(phi_nom); sin(phi_nom) cos(phi_nom)];
    e_nom_rot = a_chief * Rz_phi * vcol(e_nom_hat);
    
    % i-vector doesn't rotate in ROE space
    i_nom_rot = a_chief * vcol(i_nom_hat);
    
    % Compute errors relative to rotating/static nominals
    err_hist(idx,:) = vrow(roe_current_m) - vrow(roe_target);
    e_norm(idx) = norm(roe_current_m(3:4) - e_nom_rot);
    i_norm(idx) = norm(roe_current_m(5:6) - i_nom_rot);
end

% Also update cumulative Δv history for smooth plotting
cum_dv_hist = zeros(num_points, 1);
if ~isempty(burn_t)
    for i = 1:num_points
        burns_before = burn_t <= t_grid(i);
        if any(burns_before)
            cum_dv_hist(i) = sum(burn_dv(burns_before));
        end
    end
end

% 1) ROE-tracking error (FIXED)
figure('Name','Control-tracking error','Color','w');
ax = tiledlayout(3,2,"TileSpacing","compact");
titles = {'aδa','aδλ','aδe_x','aδe_y','aδi_x','aδi_y'};
for k = 1:6
    nexttile(ax); hold on
    plot(t_orbit, err_hist(:,k), 'LineWidth',1);
    yline(0,'k:');
    ylabel('[m]'); title(titles{k});
    grid on
    
    % Add maneuver markers
    if ~isempty(maneuver_log)
        for m = 1:length(maneuver_log)
            xline(maneuver_log(m).t_det/T, 'r:', 'Alpha', 0.3);
        end
    end
end
xlabel(ax,'Orbit');

% 2) Control window visualization (FIXED)
figure('Name','Control Window Status','Color','w');
subplot(2,1,1);
plot(t_orbit, e_norm, 'b', 'LineWidth', 1.5); hold on;
yline(delta_e_max_m, 'r--', 'LineWidth', 1.5);
fill([t_orbit(1) t_orbit(end) t_orbit(end) t_orbit(1)], ...
     [0 0 delta_e_max_m delta_e_max_m], 'r', 'FaceAlpha', 0.1);
     
% Add maneuver markers
if ~isempty(maneuver_log)
    for m = 1:length(maneuver_log)
        xline(maneuver_log(m).t_det/T, 'k:', 'Alpha', 0.5);
    end
end

ylabel('[m]'); title('Eccentricity vector error from rotating nominal');
legend('||δe - δe_{nom,rot}||', 'Control limit', 'Maneuvers', 'Location', 'best');
grid on; 
ylim([0, max(max(e_norm)*1.1, 1.5*delta_e_max_m)]);

subplot(2,1,2);
plot(t_orbit, i_norm, 'b', 'LineWidth', 1.5); hold on;
yline(delta_i_max_m, 'r--', 'LineWidth', 1.5);
fill([t_orbit(1) t_orbit(end) t_orbit(end) t_orbit(1)], ...
     [0 0 delta_i_max_m delta_i_max_m], 'r', 'FaceAlpha', 0.1);
     
% Add maneuver markers
if ~isempty(maneuver_log)
    for m = 1:length(maneuver_log)
        xline(maneuver_log(m).t_det/T, 'k:', 'Alpha', 0.5);
    end
end

ylabel('[m]'); title('Inclination vector error from static nominal');
legend('||δi - δi_{nom}||', 'Control limit', 'Maneuvers', 'Location', 'best');
grid on; 
ylim([0, max(max(i_norm)*1.1, 1.5*delta_i_max_m)]);
xlabel('Orbit');

% 3) Maneuver schedule & Δv budget (FIXED)
figure('Name','Maneuver schedule & Δv budget','Color','w');

% Left axis: Individual maneuver Δv
yyaxis left
if ~isempty(burn_t)
    stem(burn_t/T, burn_dv, 'filled', 'MarkerSize', 6);
end
ylabel('|Δv| per burn [m/s]');
ylim([0, max(burn_dv)*1.2]);

% Right axis: Cumulative Δv
yyaxis right
plot(t_orbit, cum_dv_hist, 'LineWidth', 2);
ylabel('Cumulative Δv [m/s]');
xlabel('Orbit');
grid on; 
title('Maneuver epochs and fuel use');
ylim([0, max(cum_dv_hist)*1.1]);

% 4) Δv direction statistics
if ~isempty(burn_vecs) && size(burn_vecs, 1) > 0
    dv_split = sum(abs(burn_vecs), 1);
else
    dv_split = [0, 0, 0];
end

figure('Name','Δv direction statistics','Color','w');
bar(categorical({'R', 'T', 'N'}), dv_split);
ylabel('∑|Δv| [m/s]'); grid on;
title('Where the fuel goes');

% 5) e/i Vector Phase Plot (FIXED)
figure('Name','e/i Vector Evolution','Color','w');
theta = linspace(0, 2*pi, 100);

subplot(1,2,1);
% Plot actual trajectory
plot(rel_oe(:,3), rel_oe(:,4), 'b', 'LineWidth', 1.5); hold on;

% Plot control windows at several time snapshots
time_snapshots = linspace(0, tend, 10);
cmap = cool(length(time_snapshots));
for j = 1:length(time_snapshots)
    % Recalculate phi_dot at this time (it may change slightly)
    idx_j = find(t_grid >= time_snapshots(j), 1, 'first');
    if isempty(idx_j), idx_j = num_points; end
    
    eta_j = sqrt(1 - TSX_oe(idx_j,2)^2);
    gamma_j = 0.5*J2*(Re/a_chief)^2*(1/(eta_j^4));
    phi_dot_j = 1.5*gamma_j*(5*cos(TSX_oe(idx_j,3))^2 - 1);
    phi_j = phi_dot_j * time_snapshots(j);
    
    Rz_j = [cos(phi_j) -sin(phi_j); sin(phi_j) cos(phi_j)];
    e_center_j = a_chief * Rz_j * vcol(e_nom_hat);
    
    if j == 1 || j == length(time_snapshots)
        % Full circle for first and last
        plot(e_center_j(1) + delta_e_max_m*cos(theta), ...
             e_center_j(2) + delta_e_max_m*sin(theta), '--', ...
             'Color', cmap(j,:), 'LineWidth', 1.5);
    else
        % Just a dot for intermediate positions
        plot(e_center_j(1), e_center_j(2), 'o', ...
             'Color', cmap(j,:), 'MarkerSize', 8, 'MarkerFaceColor', cmap(j,:));
    end
end

% Mark maneuver points
if ~isempty(maneuver_log)
    for k = 1:length(maneuver_log)
        idx_k = maneuver_log(k).idx;
        if idx_k <= size(rel_oe,1)
            plot(rel_oe(idx_k,3), rel_oe(idx_k,4), 'ro', ...
                 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        end
    end
end

xlabel('a·δe_x [m]'); ylabel('a·δe_y [m]');
title('Relative e-vector (control window rotates)'); 
axis equal; grid on;
legend('Trajectory', 'Initial window', 'Final window', 'Window centers', 'Maneuvers', ...
       'Location', 'best');

subplot(1,2,2);
% Plot i-vector (static control window)
plot(roe_target(5) + delta_i_max_m*cos(theta), ...
     roe_target(6) + delta_i_max_m*sin(theta), 'r--', 'LineWidth', 1.5);
hold on;

% Plot actual trajectory
plot(rel_oe(:,5), rel_oe(:,6), 'b', 'LineWidth', 1.5);

% Mark maneuver points
if ~isempty(maneuver_log)
    for k = 1:length(maneuver_log)
        idx_k = maneuver_log(k).idx;
        if idx_k <= size(rel_oe,1)
            plot(rel_oe(idx_k,5), rel_oe(idx_k,6), 'ro', ...
                 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        end
    end
end

xlabel('a·δi_x [m]'); ylabel('a·δi_y [m]');
title('Relative i-vector (static control window)'); 
axis equal; grid on;
legend('Control window', 'Trajectory', 'Maneuvers', 'Location', 'best');

% 6) ROE evolution over time
figure('Name','ROE Evolution','Color','w');
subplot(2,1,1);
plot(t_orbit, rel_oe(:,1), 'b', 'LineWidth', 1.5); hold on;
plot(t_orbit, rel_oe(:,2), 'r', 'LineWidth', 1.5);

% Add maneuver markers
if ~isempty(maneuver_log)
    for m = 1:length(maneuver_log)
        xline(maneuver_log(m).t_det/T, 'k:', 'Alpha', 0.3);
    end
end

ylabel('[m]'); title('Along-track ROE components');
legend('a·δa', 'a·δλ', 'Location', 'best'); grid on;

subplot(2,1,2);
plot(t_orbit, rel_oe(:,3), 'b', 'LineWidth', 1.5); hold on;
plot(t_orbit, rel_oe(:,4), 'r', 'LineWidth', 1.5);
plot(t_orbit, rel_oe(:,5), 'g', 'LineWidth', 1.5);
plot(t_orbit, rel_oe(:,6), 'm', 'LineWidth', 1.5);

% Add maneuver markers
if ~isempty(maneuver_log)
    for m = 1:length(maneuver_log)
        xline(maneuver_log(m).t_det/T, 'k:', 'Alpha', 0.3);
    end
end

ylabel('[m]'); title('Cross-track ROE components');
legend('a·δe_x', 'a·δe_y', 'a·δi_x', 'a·δi_y', 'Location', 'best'); grid on;
xlabel('Orbit');

% 7) RTN relative motion plots
figure('Name','RTN Relative Motion','Color','w');

% TR plane
subplot(1,3,1);
plot(rT, rR, '-', 'LineWidth', 1.5); grid on; axis equal;
xlabel('r_T [m]'); ylabel('r_R [m]');
title('RTN: TR plane');

% NR plane
subplot(1,3,2);
plot(rN, rR, '-', 'LineWidth', 1.5); grid on; axis equal;
xlabel('r_N [m]'); ylabel('r_R [m]');
title('RTN: NR plane');

% TN plane
subplot(1,3,3);
plot(rT, rN, '-', 'LineWidth', 1.5); grid on; axis equal;
xlabel('r_T [m]'); ylabel('r_N [m]');
title('RTN: TN plane');

% 8) Control Performance Analysis (NEW)
figure('Name','Control Performance Analysis','Color','w');

% Plot 1: e-vector error with control windows
subplot(2,2,1);
plot(t_orbit, e_norm, 'b', 'LineWidth', 1.5); hold on;
yline(delta_e_max_m, 'r--', 'LineWidth', 1.5);

% Mark maneuver times
if ~isempty(maneuver_log)
    for k = 1:length(maneuver_log)
        xline(maneuver_log(k).t_det/T, 'k:', 'Alpha', 0.3);
    end
end

ylabel('[m]'); xlabel('Orbit');
title('||δe - δe_{nom,rot}|| vs time');
legend('Error', 'Limit', 'Maneuvers', 'Location', 'best');
grid on; ylim([0, max(e_norm)*1.1]);

% Plot 2: i-vector error with control windows
subplot(2,2,2);
plot(t_orbit, i_norm, 'b', 'LineWidth', 1.5); hold on;
yline(delta_i_max_m, 'r--', 'LineWidth', 1.5);

% Mark maneuver times
if ~isempty(maneuver_log)
    for k = 1:length(maneuver_log)
        xline(maneuver_log(k).t_det/T, 'k:', 'Alpha', 0.3);
    end
end

ylabel('[m]'); xlabel('Orbit');
title('||δi - δi_{nom}|| vs time');
legend('Error', 'Limit', 'Maneuvers', 'Location', 'best');
grid on; ylim([0, max(i_norm)*1.1]);

% Plot 3: Maneuver breakdown
subplot(2,2,3);
if ~isempty(maneuver_log)
    orbits_man = [maneuver_log.t_det]/T;
    dvt1_array = abs([maneuver_log.dvt1]);
    dvt2_array = abs([maneuver_log.dvt2]);
    dvn_array = [maneuver_log.dvn];
    
    bar(orbits_man, [dvt1_array; dvt2_array; dvn_array]', 'stacked');
    xlabel('Orbit'); ylabel('Δv [m/s]');
    title('Maneuver breakdown');
    legend('|Δv_{T1}|', '|Δv_{T2}|', 'Δv_N', 'Location', 'best');
    grid on;
    xlim([0, tend/T]);
end

% Plot 4: Error phase space
subplot(2,2,4);
plot(e_norm, i_norm, 'b-', 'LineWidth', 1); hold on;
plot(e_norm(1), i_norm(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(e_norm(end), i_norm(end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

% Draw control box
rectangle('Position', [0, 0, delta_e_max_m, delta_i_max_m], ...
          'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');

xlabel('||δe - δe_{nom,rot}|| [m]'); 
ylabel('||δi - δi_{nom}|| [m]');
title('Error phase space');
legend('Trajectory', 'Start', 'End', 'Control window', 'Location', 'best');
grid on; axis equal;
xlim([0, max(e_norm)*1.1]); ylim([0, max(i_norm)*1.1]);

% Fixed Error Evolution Analysis section
fprintf('\n=== ERROR EVOLUTION ANALYSIS (FIXED) ===\n');
fprintf('%-12s %-12s %-12s %-12s %-12s\n', 'Orbit', 'de_error', 'di_error', '||de||', '||di||');
fprintf('%-12s %-12s %-12s %-12s %-12s\n', '-----', '--------', '--------', '------', '------');

key_orbits = [0, 1, 5, 10, 15, 20, 25, 30, 40, 45];
for orbit_num = key_orbits
    idx = find(t_orbit >= orbit_num, 1, 'first');
    if ~isempty(idx) && idx <= num_points
        % Use stored e_norm and i_norm directly - they already contain the errors
        fprintf('%-12.1f %-12.2f %-12.2f %-12.2f %-12.2f\n', ...
                orbit_num, e_norm(idx), i_norm(idx), ...
                norm(rel_oe(idx,3:4)), norm(rel_oe(idx,5:6)));
    end
end

% Print final summary with more detail
fprintf('\n=== FINAL SUMMARY ===\n');
fprintf('Simulation duration: %.1f days (%.1f orbits)\n', sim_days, tend/T);
fprintf('Total maneuvers: %d\n', length(maneuver_log));
fprintf('Total Δv: %.3f m/s\n', cum_dv);
fprintf('Average Δv per maneuver: %.4f m/s\n', cum_dv/max(1, length(maneuver_log)));
fprintf('Average Δv per day: %.4f m/s\n', cum_dv/sim_days);

if ~isempty(burn_vecs)
    fprintf('\nΔv breakdown by direction:\n');
    fprintf('  Radial (R):     %.4f m/s (%.1f%%)\n', dv_split(1), 100*dv_split(1)/sum(dv_split));
    fprintf('  Tangential (T): %.4f m/s (%.1f%%)\n', dv_split(2), 100*dv_split(2)/sum(dv_split));
    fprintf('  Normal (N):     %.4f m/s (%.1f%%)\n', dv_split(3), 100*dv_split(3)/sum(dv_split));
end

fprintf('\nFinal errors:\n');
fprintf('  δe error: %.2f m (%.1f%% of limit)\n', e_norm(end), 100*e_norm(end)/delta_e_max_m);
fprintf('  δi error: %.2f m (%.1f%% of limit)\n', i_norm(end), 100*i_norm(end)/delta_i_max_m);

%% Helper Functions

%% Fixed Helper Functions

% Fixed Helper Functions

function dirECI = rtn_axis(oeChief, oeDep, which)
    % Return unit R / T / N direction (in ECI) at the CHIEF position
    % Returns a COLUMN vector for consistency
    global mu tol
    
    % Get CHIEF position and velocity
    rvC = oe2rv([oeChief(1:5)  mean2true(oeChief(6), oeChief(2), tol)], mu);
    
    % Use CHIEF state to define RTN frame
    [R,T,N] = eci2rtn_dir(rvC(1:3), rvC(4:6));
    
    switch upper(which)
        case 'R', dirECI = R(:);  % Force column vector
        case 'T', dirECI = T(:);  % Force column vector
        otherwise , dirECI = N(:); % Force column vector
    end
end

function [R,T,N] = eci2rtn_dir(r, v)
    % Returns row vectors
    R = r(:).' / norm(r);           % radial (row)
    h = cross(r, v);
    N = h(:).' / norm(h);           % normal (row)
    T = cross(N, R);                % transverse (row)
end

function oeNew = kick(oeRow, dvECI)
    % Add an impulsive Δv (ECI column vector) to an OE row; return OE row
    global mu tol
    rv = oe2rv([oeRow(1:5) mean2true(oeRow(6), oeRow(2), tol)], mu);
    rv(4:6) = rv(4:6) + dvECI(:);  % Ensure column vector addition
    oe = safe_rv2oe(rv, mu);
    oeNew = [oe(1:5) safe_true2mean(oe(6), oe(2))];
end