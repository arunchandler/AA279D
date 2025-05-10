clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

%% Problem 1

scenario = input( ...
    ['Select rel_qns_init scenario:\n' ...
     '  1: DEM formation (e_term=122/√2)\n' ...
     '  2: Example “circular” formation\n' ...
     '  3: Custom (enter manually)\n' ...
     'Your choice: '] );

switch scenario
    case 1
        % DEM
        e_term = 122/sqrt(2);
        rel_qns_init = [0, 340, e_term, e_term, 0, 256];
        scenario_name = 'DEM';
        
    case 2
        % Pursuit
        e_term = 4982/sqrt(2);
        rel_qns_init = [0, -76050, e_term, e_term, 80, 0 ];  
        scenario_name = 'Pursuit';
        
    case 3
        % Large cross-track
        e_term = 3600/sqrt(2);
        rel_qns_init = [0, 0, e_term, e_term, 250, 0 ];  
        scenario_name = 'Large cross-track';
    
    case 4
        % Short baseline
        e_term = 250/sqrt(2);
        rel_qns_init = [0, -340, e_term, e_term, 250, 0 ];  
        scenario_name = 'Short baseline';

    case 5
        %Damico separation
        rel_qns_init = [0, 0, 0, 300, 0, -1000];
        scenario_name = 'Damico Paper Initialization';

        
    otherwise
        error('Invalid scenario');
end

fprintf('Running scenario "%s": rel_qns_init = [%g %g %g %g %g %g]\n\n', ...
        scenario_name, rel_qns_init);

% initial chief elements
a_TSX_init    = 6892944.65;  % m
e_TSX_init    = 0.00013700;
i_TSX_init    = deg2rad(97.440124);
RAAN_TSX_init = deg2rad(104.274891);
omega_TSX_init= deg2rad(67.975723);
M_TSX_init    = deg2rad(292.169756);
nu_TSX_init   = mean2true(M_TSX_init, e_TSX_init, tol);
u_TSX_init    = nu_TSX_init + omega_TSX_init;

TSX_init_oe = [a_TSX_init, e_TSX_init, i_TSX_init, ...
               RAAN_TSX_init, omega_TSX_init, M_TSX_init];

% map to TDX via your chosen rel_qns_init
TDX_init_oe = qns2oe(TSX_init_oe, rel_qns_init);

% timing parameters, propagation, etc. (unchanged)
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init^3);
T          = 2*pi/n;
n_orbit = ceil((30 * 86400) / T);  % Number of orbits in 30 days
tend = 30 * 86400;  % 30 days in seconds
num_points = 10000;  % Increase for better resolution
dt = (tend - tstart)/(num_points-1);
t_grid = linspace(tstart, tend, num_points).';
t_orbit = t_grid / T;

% propagate with your ode4 + J2
[t_out, TSX_oe] = ode4(@compute_rates_GVE_J2, [tstart, tend]', TSX_init_oe', dt);
[~,     TDX_oe] = ode4(@compute_rates_GVE_J2, [tstart, tend]', TDX_init_oe,  dt);

TSX_oe(:,6) = wrapTo2Pi(TSX_oe(:,6));
TDX_oe(:,6) = wrapTo2Pi(TDX_oe(:,6));

% compute relative OEs and RTN
rel_oe = zeros(num_points, 6);
rtn    = zeros(num_points, 6);

for idx = 1:num_points
    
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
end

%plotting orbital elements
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

% unpack
rR = rtn(:,1);
rT = rtn(:,2);
rN = rtn(:,3);

% TR plane: T vs R
% NR plane: N vs R
% TN plane: T vs N
figure;
subplot(1,3,1)
plot(rT, rR, '-'); grid on; axis equal
xlabel('r_T [m]'); ylabel('r_R [m]');
title('RTN: TR plane')

subplot(1,3,2)
plot(rN, rR, '-'); grid on; axis equal
xlabel('r_N [m]'); ylabel('r_R [m]');
title('RTN: NR plane')

subplot(1,3,3)
plot(rT, rN, '-'); grid on; axis equal
xlabel('r_T [m]'); ylabel('r_N [m]');
title('RTN: TN plane')

%Now using J2 STM matrix
roe_stm = zeros(num_points, 6);
roe_stm(1,:) = rel_qns_init;

for idx = 2:num_points
 
    oe = TSX_oe(idx,:);
    phi = stm_qns_j2(dt, oe);
    roe_stm(idx,:) = phi * roe_stm(idx-1,:)';

end

% Unpack propagated relative OE
delta_a_prop      = rel_oe(:,1);
delta_lambda_prop = rel_oe(:,2);
delta_ex_prop     = rel_oe(:,3);
delta_ey_prop     = rel_oe(:,4);
delta_ix_prop     = rel_oe(:,5);
delta_iy_prop     = rel_oe(:,6);

% Unpack STM relative OE
delta_a_stm      = roe_stm(:,1);
delta_lambda_stm = roe_stm(:,2);
delta_ex_stm     = roe_stm(:,3);
delta_ey_stm     = roe_stm(:,4);
delta_ix_stm     = roe_stm(:,5);
delta_iy_stm     = roe_stm(:,6);

% 1) Δa vs Δλ
figure;
plot(delta_a_prop,      delta_lambda_prop, '-' , 'LineWidth',1.5); hold on;
plot(delta_a_stm,       delta_lambda_stm, '--', 'LineWidth',1.5);
grid on; axis equal;
xlabel('a\delta a [m]'); ylabel('a\delta \lambda [m]');
title('Relative OE: a\delta a vs a\delta \lambda');
legend('Propagated','STM','Location','best');

% 2) Δe_x vs Δe_y
figure;
plot(delta_ex_prop,     delta_ey_prop, '-' , 'LineWidth',1.5); hold on;
plot(delta_ex_stm,      delta_ey_stm, '--', 'LineWidth',1.5);
grid on; axis equal;
xlabel('a\delta e_x [m]'); ylabel('a\delta e_y [m]');
title('Relative OE: a\delta e_x vs a\delta e_y');
legend('Propagated','STM','Location','best');

% 3) Δi_x vs Δi_y
figure;
plot(delta_ix_prop,     delta_iy_prop, '-' , 'LineWidth',1.5); hold on;
plot(delta_ix_stm,      delta_iy_stm, '--', 'LineWidth',1.5);
grid on; axis equal;
xlabel('a\delta i_x [m]'); ylabel('a\delta i_y [m]');
title('Relative OE: a\delta i_x vs a\delta i_y');
legend('Propagated','STM','Location','best');


% === CONTROL SECTION: Initialize ===
a = a_TSX_init;
v = sqrt(mu / a);
maneuver_log = [];

% Target (nominal) ROE [δa, δλ, δeₓ, δeᵧ, δiₓ, δiᵧ] in meters
roe_target = [0; 0; 0; 300; 0; -1000];

% Maneuver detection threshold (angle between e/i vectors)
angle_thresh_deg = 7;  % degrees

% === Loop through to detect maneuver triggers and compute Δv ===
for idx = 1:num_points
    roe_now = rel_oe(idx,:)';
    delta_e_vec = roe_now(3:4);
    delta_i_vec = roe_now(5:6);

    % Angle between e and i vectors
    cos_angle = dot(delta_e_vec, delta_i_vec) / (norm(delta_e_vec) * norm(delta_i_vec));
    angle_rad = acos(max(min(cos_angle, 1), -1));  % Clamp for numerical safety
    angle_deg = rad2deg(angle_rad);
    sep_deg = min(angle_deg, abs(180 - angle_deg));  % True angular separation    
    if mod(t_orbit(idx), 10) < dt/T  % approximately once per orbit
        fprintf('[%.2f orbits] e/i angle separation = %.2f°\n', t_orbit(idx), sep_deg);
    end

    % Check if vectors are sufficiently misaligned AND no prior trigger
    if sep_deg > angle_thresh_deg && isempty(maneuver_log)
        % Normalize error terms (unitless)
        delta_e_err = (delta_e_vec - roe_target(3:4)) / a;
        delta_a_err = (roe_now(1) - roe_target(1)) / a;

        delta_e_mag = norm(delta_e_err);

        % Δv computation
        dvt1 = (v/4)*(delta_e_mag + delta_a_err);
        dvt2 = -(v/4)*(delta_e_mag - delta_a_err);

        % Maneuver phasing (location)
        u1 = atan2(delta_e_err(2), delta_e_err(1));
        u2 = mod(u1 + pi, 2*pi);

        % Log maneuver info
        maneuver_log = [maneuver_log; struct( ...
            'idx', idx, ...
            't', t_grid(idx), ...
            'dvt1', dvt1, ...
            'dvt2', dvt2, ...
            'u1', u1, ...
            'u2', u2, ...
            'sep_deg', sep_deg, ...
            'delta_e_err', delta_e_err)];

        fprintf('[%.2f orbits] Trigger: Δv_T1 = %.4f m/s, Δv_T2 = %.4f m/s, separation = %.2f°\n', ...
            t_orbit(idx), dvt1, dvt2, sep_deg);

        
    end
end

% === APPLY TWO-PULSE MANEUVER TO TDX ===
if ~isempty(maneuver_log)
    t1   = maneuver_log(1).t;
    dvt1 = maneuver_log(1).dvt1;
    dvt2 = maneuver_log(1).dvt2;

    % re‐compute half‐period
    n       = sqrt(mu / a_TSX_init^3);
    T       = 2*pi / n;
    dt_half = T / 2;

    % === First Impulse ===
    [~, idx1] = min(abs(t_grid - t1));
    oe1_mean  = TDX_oe(idx1, :);            % [a, e, i, Ω, ω, M_mean]

    % unpack and convert M→ν
    a1      = oe1_mean(1);
    e1      = oe1_mean(2);
    i1      = oe1_mean(3);
    RAAN1   = oe1_mean(4);
    omega1  = oe1_mean(5);
    M1_mean = oe1_mean(6);
    nu1     = mean2true(M1_mean, e1, tol);

    % rebuild with true anomaly
    oe1_true = [a1, e1, i1, RAAN1, omega1, nu1];

    % get r1, v1
    rv1 = oe2rv(oe1_true, mu);
    r1  = rv1(1:3);
    v1  = rv1(4:6);

    % RTN frame at impulse 1
    h1     = cross(r1, v1);
    R_hat  = r1 / norm(r1);
    N_hat  = h1 / norm(h1);
    T_hat  = cross(N_hat, R_hat);
    dv1_vec = dvt1 * T_hat;

    % apply Δv1
    v1_new = v1 + dv1_vec;

    % back to OE (rv→oe yields ν, so convert back to M)
    X1      = [r1; v1_new];
    oe1_full = rv2oe(X1, mu);      % returns [a,e,i,Ω,ω,ν]
    nu1_new = oe1_full(6);
    e1_new  = oe1_full(2);
    M1_new  = true2mean(nu1_new, e1_new);
    oe1_new = [oe1_full(1:5), M1_new];

    % === Propagate to Second Impulse ===
    t2 = t1 + dt_half;
    [~, idx2_start] = min(abs(t_grid - t1));
    [~, idx2_end]   = min(abs(t_grid - t2));
    dt_fine = (t2 - t1) / (idx2_end - idx2_start + 1);

    [~, oe_between] = ode4(@compute_rates_GVE_J2, [t1 t2]', oe1_new', dt_fine);

    % === Second Impulse ===
    oe2_mean = oe_between(end, :);        % [a, e, i, Ω, ω, M_mean]

    % unpack and convert M→ν
    a2       = oe2_mean(1);
    e2       = oe2_mean(2);
    i2       = oe2_mean(3);
    RAAN2    = oe2_mean(4);
    omega2   = oe2_mean(5);
    M2_mean  = oe2_mean(6);
    nu2      = mean2true(M2_mean, e2, tol);

    oe2_true = [a2, e2, i2, RAAN2, omega2, nu2];

    % get r2, v2
    rv2 = oe2rv(oe2_true, mu);
    r2  = rv2(1:3);
    v2  = rv2(4:6);

    % RTN at impulse 2
    h2      = cross(r2, v2);
    R_hat2  = r2 / norm(r2);
    N_hat2  = h2 / norm(h2);
    T_hat2  = cross(N_hat2, R_hat2);
    dv2_vec = dvt2 * T_hat2;

    % apply Δv2
    v2_new = v2 + dv2_vec;

    % back to OE
    X2       = [r2; v2_new];
    oe2_full = rv2oe(X2, mu);       % [a,e,i,Ω,ω,ν]
    nu2_new = oe2_full(6);
    e2_new  = oe2_full(2);
    M2_new  = true2mean(nu2_new, e2_new);
    oe2_new = [oe2_full(1:5), M2_new];

    % === Final Propagation to End ===
    [~, oe_after] = ode4(@compute_rates_GVE_J2, [t2 tend]', oe2_new', dt);

    % === Stitch Everything Back Together ===
    TDX_oe = [ TDX_oe(1:idx1-1, :);
               oe1_new;
               oe_between(2:end, :);
               oe2_new;
               oe_after(2:end, :) ];

    % pick out post‑burn rel‑OE (in meters)
    rel_post = TSX_oe(idx2_end,1) * compute_roes( ...
                TSX_oe(idx2_end,:), ...
                TDX_oe(idx2_end,:) )';

    % split out e‑ and i‑vectors
    delta_e = rel_post(3:4);
    delta_i = rel_post(5:6);

    % magnitudes
    norm_e = norm(delta_e);
    norm_i = norm(delta_i);

    % e/i separation angle
    cos_sep = dot(delta_e, delta_i) / (norm_e * norm_i);
    sep_deg_post = rad2deg(acos( max(min(cos_sep,1),-1) ));

    fprintf('\nPost‑burn summary:\n');
    fprintf('   ||δe|| = %8.2f m   (target %8.2f m)\n', norm_e, norm(roe_target(3:4)));
    fprintf('   ||δi|| = %8.2f m   (target %8.2f m)\n', norm_i, norm(roe_target(5:6)));
    fprintf('   e/i separation = %6.2f°   (target 180°)\n\n', sep_deg_post);

    % report total Δv
    total_dv = norm(dv1_vec) + norm(dv2_vec);
    fprintf('\n>>> Total Δv applied: %.6f m/s (Δv₁ + Δv₂)\n\n', total_dv);
end