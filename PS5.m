clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

% Problem 1 - relative dynamics without control

scenario = input( ...
    ['Select rel_qns_init scenario:\n' ...
     '  1: M-D2\n' ...
     '  2: M-D3\n' ...
     '  3: M-D4\n' ...
     'Your choice: '] );

switch scenario
    case 1
        % M-D2
        rel_qns_init = [0, 0, 0, 300, 0, 400];
        scenario_name = 'M-C';
        
    case 2
        % M-D3
        rel_qns_init = [0, 0, 0, 300, 0, 500];  
        scenario_name = 'M-D1';
        
    case 3
        % M-D4
        rel_qns_init = [0, 0, 0, 500, 0, 300];  
        scenario_name = 'M-E';
        
    otherwise
        error('Invalid scenario');
end

fprintf('Running scenario "%s": rel_qns_init = [%g %g %g %g %g %g]\n\n', ...
        scenario_name, rel_qns_init);

% initial chief elements
a_TSX_init    = 6886536.686;  % m
e_TSX_init    = 0.0001264;
i_TSX_init    = deg2rad(97.4453);
RAAN_TSX_init = deg2rad(351.0108);
omega_TSX_init= deg2rad(101.2452);
M_TSX_init    = deg2rad(11.6520);
nu_TSX_init   = mean2true(M_TSX_init, e_TSX_init, tol);
u_TSX_init    = M_TSX_init + omega_TSX_init;

TSX_init_oe = [a_TSX_init, e_TSX_init, i_TSX_init, ...
               RAAN_TSX_init, omega_TSX_init, M_TSX_init];

% map to TDX via your chosen rel_qns_init
TDX_init_oe = qns2oe(TSX_init_oe, rel_qns_init);

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

% Problem 2 - reconfiguration

scenario = input( ...
    ['Select reconfiguration scenario:\n' ...
     '  1: M-D2 to M-D3\n' ...
     '  2: M-D3 to M-D4\n' ...
     '  3: Test Reconfiguration\n' ...
     '  4: Test Formation Keeping' ...
     'Your choice: '] );

switch scenario
    case 1
        % M-D2 to M-D3
        rel_qns_pre = [0, 0, 0, 300, 0, 400];
        rel_qns_post  = [0, 0, 0, 300, 0, 500];
        scenario_name = 'M-D2 to M-D3';
        
    case 2
        % M-D3 to M-D4
        rel_qns_pre = [0, 0, 0, 300, 0, 500];
        rel_qns_post  = [0, 0, 0, 500, 0, 300];
        scenario_name = 'M-D3 to M-D4';

    case 3
        % Test
        rel_qns_pre = [0, 0, 500, 500, 500, 500];
        rel_qns_post  = [0, 0, 300, 300, 300, 300];
        scenario_name = 'Test Reconfiguration';

    case 4
        % Test Formation Keeping
        rel_qns_pre = [0, 0, 500, 500, 500, 500];
        rel_qns_post = [0, 0, 500, 500, 500, 500];
        senario_name = 'Test Formation Keeping';
        
    otherwise
        error('Invalid scenario');
end

fprintf('Running scenario "%s": rel_qns_init = [%g %g %g %g %g %g]\n rel_qns_fin = [%g %g %g %g %g %g]\n', ...
        scenario_name, [rel_qns_pre, rel_qns_post]);

TDX_init_oe_2 = qns2oe(TSX_init_oe, rel_qns_pre);
TSX_init_rv = oe2rv(TSX_init_oe, mu);
TDX_init_rv = oe2rv(TDX_init_oe_2, mu);
TDX_init_rtn = eci2rtn(TSX_init_rv, TDX_init_rv);
rel_state_init = [TDX_init_rtn; TSX_init_rv];

state_out = zeros(num_points, 12);
state_out(1,:) = rel_state_init';

TSX_oe_2 = zeros(num_points, 6);
TDX_oe_2 = zeros(num_points, 6);
TSX_oe_2(1,:) = TSX_init_oe;
TDX_oe_2(1,:) = TDX_init_oe_2;

rel_oe_2 = zeros(num_points, 6);
rel_oe_2(1,:) = a_TSX_init*compute_roes(TSX_init_oe, TDX_init_oe_2);

state_cur = rel_state_init;

% ----naive least squares----

%delta v time calculations
n_man = 3;

delta_v_times = linspace(tstart, tend, n_man+2);
delta_v_times = delta_v_times(2:end-1)
delta_v_vals = naive_least_squares(delta_v_times, rel_qns_pre./a_TSX_init, rel_qns_post./a_TSX_init, TSX_init_oe, u_TSX_init, tstart, tend)

mask = (orbit_num == orbit_num);

DV1 = norm(delta_v_vals(1,:));
DV2 = norm(delta_v_vals(2,:));
DV3 = norm(delta_v_vals(3,:));
DVT = DV1 + DV2 + DV3;
DV = sqrt( sum(delta_v_vals.^2, 2) );

tmp = find(mask);  % list of indices
[~,loc]    = min(abs(t_grid - delta_v_times(1)));
idx_1    = tmp(loc);

[~,loc]    = min(abs(t_grid - delta_v_times(2)));
idx_2    = tmp(loc);

[~,loc]    = min(abs(t_grid - delta_v_times(3)));
idx_3     = tmp(loc);

for idx = 2:num_points

    t_cur = t_grid(idx-1);
    t_next = t_grid(idx);

    TSX_ECI_cur = state_cur(7:12);
    TDX_RTN_cur = state_cur(1:6);

    % apply ΔV’s at the target indices
    if idx == idx_1
        TDX_RTN_cur(4:6) = TDX_RTN_cur(4:6) + delta_v_vals(1,:)';
    end
    if idx == idx_2
        TDX_RTN_cur(4:6) = TDX_RTN_cur(4:6) + delta_v_vals(2,:)';
    end
    if idx == idx_3
        TDX_RTN_cur(4:6) = TDX_RTN_cur(4:6) + delta_v_vals(3,:)';
    end

    % write the modified RTN back into state_cur
    [TDX_ECI_cur, ~] = rtn2eci(TSX_ECI_cur, TDX_RTN_cur);

    [~, TDX_ECI_next] = ode4(@compute_rates_rv_perturbed, [t_cur, t_next]', TDX_ECI_cur,  dt);
    [~, TSX_ECI_next] = ode4(@compute_rates_rv_perturbed, [t_cur, t_next]', TSX_ECI_cur,  dt);
    TDX_ECI_next = TDX_ECI_next(2,:)';
    TSX_ECI_next = TSX_ECI_next(2,:)';
    [TDX_RTN_next,~] = eci2rtn(TSX_ECI_next, TDX_ECI_next);
    state_next = [TDX_RTN_next; TSX_ECI_next];
    state_out(idx,:) = state_next;

    TSX_params = rv2oe(TSX_ECI_next, mu);
    TSX_oe_2(idx,:) = [TSX_params(1:5), true2mean(TSX_params(6), TSX_params(2))];
    
    a1 = TSX_oe_2(idx,1); e1 = TSX_oe_2(idx,2); i1 = TSX_oe_2(idx,3);
    RAAN1 = TSX_oe_2(idx,4); omega1 = TSX_oe_2(idx,5);
    M1 = wrapTo2Pi(TSX_oe_2(idx,6)); u1 = M1 + omega1;

    TDX_params = rv2oe(TDX_ECI_next, mu);
    TDX_oe_2(idx,:) = [TDX_params(1:5), true2mean(TDX_params(6), TDX_params(2))];
    
    a2 = TDX_oe_2(idx,1); e2 = TDX_oe_2(idx,2); i2 = TDX_oe_2(idx,3);
    RAAN2 = TDX_oe_2(idx,4); omega2 = TDX_oe_2(idx,5);
    M2 = wrapTo2Pi(TDX_oe_2(idx,6)); u2 = M2 + omega2;
    
    rel_oe_2(idx,:) = a1 * compute_roes([a1, e1, i1, RAAN1, omega1, M1], ...
                                        [a2, e2, i2, RAAN2, omega2, M2])';

    state_cur = state_next;

end

daf = rel_oe_2(end,1);
dlambdaf = rel_oe_2(end,2);
dexf = rel_oe_2(end,3);
deyf = rel_oe_2(end,4);
dixf = rel_oe_2(end,5);
diyf = rel_oe_2(end,6);

% --- print Delta‑V budget ---
fprintf('\n=== ΔV Budget ===\n');
fprintf('DV1 = %.6f m/s\n', DV1);
fprintf('DV2 = %.6f m/s\n', DV2);
fprintf('DV3 = %.6f m/s\n', DV3);
fprintf('Total ΔV = %.6f m/s\n\n', DVT);

% --- print achieved QNS changes ---
fprintf('=== QNS Element Changes ===\n');
fprintf('%-12s  %-12s  %-12s\n','Element','Desired','Achieved');
fprintf('%-12s  % .6f   % .6f\n','Δa',rel_qns_post(1),daf);
fprintf('%-12s  % .6f   % .6f\n','Δλ',rel_qns_post(2),dlambdaf);
fprintf('%-12s  % .6f   % .6f\n','Δe_x',rel_qns_post(3),dexf);
fprintf('%-12s  % .6f   % .6f\n','Δe_y',rel_qns_post(4),deyf);
fprintf('%-12s  % .6f   % .6f\n','Δi_x',rel_qns_post(5),dixf);
fprintf('%-12s  % .6f   % .6f\n','Δi_y',rel_qns_post(6),diyf);

% ——— Visualization ———

% 1) Define split points for three segments
idxs = unique([1, idx_1, idx_2, idx_3, num_points]);
nSeg = numel(idxs) - 1;

% 2) Colormap & legend labels
colors = lines(nSeg);
legend_labels = { ...
    'Before Maneuvers', ...
    'After \Delta V1', ...
    'After \Delta V2', ...
    'After \Delta V3'
};

% ——— RTN‑frame relative motion ———
figure;

% TR plane: T on x, R on y
subplot(1,3,1); hold on;
for s = 1:nSeg
    segStart = idxs(s);
    segEnd   = min(idxs(s+1), num_points);
    if segStart < segEnd
        seg = segStart:segEnd;
        plot(state_out(seg,2), state_out(seg,1), 'Color', colors(s,:));
    end
end
hold off;
axis equal; grid on;
xlabel('T [m]'); ylabel('R [m]');
title('RTN: TR plane');
% only one legend for the entire RTN figure
legend(legend_labels, 'Location', 'best');

% NR plane: N on x, R on y
subplot(1,3,2); hold on;
for s = 1:nSeg
    seg = idxs(s):min(idxs(s+1), num_points);
    plot(state_out(seg,3), state_out(seg,1), 'Color', colors(s,:));
end
hold off;
axis equal; grid on;
xlabel('N [m]'); ylabel('R [m]');
title('RTN: NR plane');

% TN plane: T on x, N on y
subplot(1,3,3); hold on;
for s = 1:nSeg
    seg = idxs(s):min(idxs(s+1), num_points);
    plot(state_out(seg,2), state_out(seg,3), 'Color', colors(s,:));
end
hold off;
axis equal; grid on;
xlabel('T [m]'); ylabel('N [m]');
title('RTN: TN plane');

% ——— Relative QNS elements vs. Orbit number ———
figure;
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

labels = {'a\delta a [m]', ...
          'a\delta\lambda [m]', ...
          'a\delta e_x [m]', ...
          'a\delta e_y [m]', ...
          'a\delta i_x [m]', ...
          'a\delta i_y [m]'};

for k = 1:6
    nexttile;
    plot(t_orbit, rel_oe_2(:,k));
    grid on;
    xlabel('Orbit Number');
    ylabel(labels{k});
    title(labels{k});
end

% cumulative delta-v
cum_dv = zeros(size(t_orbit));
cum_dv(idx_1:end) = DV1;
cum_dv(idx_2:end) = DV1 + DV2;   
cum_dv(idx_3:end) = DV1 + DV2 + DV3;
figure;
plot(t_orbit, cum_dv);
xlabel('Orbit Number'); ylabel('Cumulative Delta-v (m/s)');
title('Cumulative Delta-v');

% --- 2D relative‐element scatter for control visualization ---
figure;
% Δλ vs Δa
subplot(1,3,1);
plot(rel_oe_2(:,2), rel_oe_2(:,1), 'LineWidth',1.2);
xlabel('a\delta\lambda [m]');
ylabel('a\deltaa [m]');
grid on;
axis equal;
% Δe_x vs Δe_y
subplot(1,3,2);
plot(rel_oe_2(:,3), rel_oe_2(:,4), 'LineWidth',1.2);
xlabel('a\deltae_x [m]');
ylabel('a\deltae_y [m]');
grid on;
axis equal;
% Δi_x vs Δi_y
subplot(1,3,3);
plot(rel_oe_2(:,5), rel_oe_2(:,6), 'LineWidth',1.2);
xlabel('a\deltai_x [m]');
ylabel('a\deltai_y [m]');
grid on;
axis equal;