clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

%% Problem 1 - relative dynamics without control

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

%% Problem 2 - reconfiguration


scenario = input( ...
    ['Select reconfiguration scenario:\n' ...
     '  1: M-D2 to M-D3\n' ...
     '  2: M-D3 to M-D4\n' ...
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

TDX_oe_2 = zeros(num_points, 6);
TDX_oe_2(1,:) = TDX_init_oe_2;

rel_oe_2 = zeros(num_points, 6);
rel_oe_2(1,:) = a_TSX_init*compute_roes(TSX_init_oe, TDX_init_oe_2);

state_cur = rel_state_init;

% --- precompute target indices ---
% chief argument of latitude array:
u_TSX = wrapTo2Pi(TSX_oe(:,6)) + TSX_oe(:,5);
target_orbit = 7;
mask = (orbit_num == target_orbit);

% delta V calculation

Ddiy     = rel_qns_post(6) - rel_qns_pre(6);
Ddey     = rel_qns_post(4) - rel_qns_pre(4);
Ddlambda = rel_qns_post(2) - rel_qns_pre(2);

a_m = TSX_oe(num_points/2,1);

%out of plane portion
uM_op = pi/2
dvn = Ddiy * n / sin(uM_op);

%in plane portion
dvt = 0.0; %because we dont want sma change, we can set this to 0 and choose dvr to satisfy change in dey
uM_ip1 = 0.0
uM_ip2 = pi - uM_ip1
dvr1 = Ddey * n * a_m / (-2 * cos(uM_ip1));
dvr2 = -dvr1;

% --- compute total delta‐V ---
DV_op  = abs(dvn);        % out‐of‐plane impulse magnitude
DV_ip1 = abs(dvr1);       % in‐plane impulse at u = 0
DV_ip2 = abs(dvr2);       % in‐plane impulse at u = pi
DV_tot = DV_op + DV_ip1 + DV_ip2;

tmp = find(mask);  % list of indices in orbit 7
[~,loc]    = min(abs(u_TSX(mask) - uM_ip1));
idx_ip1    = tmp(loc);

[~,loc]    = min(abs(u_TSX(mask) - uM_ip2));
idx_ip2    = tmp(loc);

[~,loc]    = min(abs(u_TSX(mask) - uM_op ));
idx_op     = tmp(loc);

for idx = 2:num_points

    t_cur = t_grid(idx-1);
    t_next = t_grid(idx);

    TSX_ECI_cur = state_cur(7:12);
    TDX_RTN_cur = state_cur(1:6);

    % apply ΔV’s at the target indices
    if idx == idx_ip1
        TDX_RTN_cur(4) = TDX_RTN_cur(4) + dvr1;
        fprintf('ip1\n');
    end
    if idx == idx_ip2
        TDX_RTN_cur(4) = TDX_RTN_cur(4) + dvr2;
        fprintf('ip2\n');
    end
    if idx == idx_op
        TDX_RTN_cur(6) = TDX_RTN_cur(6) + dvn;
        fprintf('op\n');
    end

    % write the modified RTN back into state_cur
    state_cur(1:6)   = TDX_RTN_cur;

    [TDX_ECI_cur, ~] = rtn2eci(TSX_ECI_cur, TDX_RTN_cur);

    [~, state_next] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [t_cur, t_next]', state_cur,  dt);
    state_next = state_next(2,:)';
    state_out(idx,:) = state_next;

    TSX_ECI_next = state_next(7:12);
    TDX_RTN_next = state_next(1:6);
    [TDX_ECI_next, R] = rtn2eci(TSX_ECI_next, TDX_RTN_next);
    
    a1 = TSX_oe(idx,1); e1 = TSX_oe(idx,2); i1 = TSX_oe(idx,3);
    RAAN1 = TSX_oe(idx,4); omega1 = TSX_oe(idx,5);
    M1 = wrapTo2Pi(TSX_oe(idx,6)); u1 = M1 + omega1;

    TDX_params = rv2oe(TDX_ECI_next, mu);
    TDX_oe_2(idx,:) = [TDX_params(1:5), true2mean(TDX_params(6), TDX_params(2))];
    
    a2 = TDX_oe_2(idx,1); e2 = TDX_oe_2(idx,2); i2 = TDX_oe_2(idx,3);
    RAAN2 = TDX_oe_2(idx,4); omega2 = TDX_oe_2(idx,5);
    M2 = wrapTo2Pi(TDX_oe_2(idx,6)); u2 = M2 + omega2;
    
    rel_oe_2(idx,:) = a1 * compute_roes([a1, e1, i1, RAAN1, omega1, M1], ...
                                        [a2, e2, i2, RAAN2, omega2, M2])';

    state_cur = state_next;

end

% --- print Delta‑V budget ---
fprintf('\n=== ΔV Budget ===\n');
fprintf('Out‑of‑plane (dvn)     = %.6f m/s\n', dvn);
fprintf('In‑plane stage 1 (dvr1) = %.6f m/s\n', dvr1);
fprintf('In‑plane stage 2 (dvr2) = %.6f m/s\n', dvr2);
fprintf('Total ΔV               = %.6f m/s\n\n', DV_tot);

% --- print achieved QNS changes ---
fprintf('=== QNS Element Changes ===\n');
fprintf('Δi_y      = %.6e rad\n', Ddiy);
fprintf('Δe_y      = %.6e    \n', Ddey);
fprintf('Δλ (dlambda) = %.6e rad\n', Ddlambda);

% ——— Visualization ———

% 1) Define split points for three segments
idxs = unique([1, idx_ip1, idx_ip2, num_points]);
nSeg = numel(idxs) - 1;

% 2) Colormap & legend labels
colors = lines(nSeg);
legend_labels = { ...
    'Before Maneuvers', ...
    'After \delta v_r_1', ...
    'After \delta v_r_2 and \delta v_t' ...
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
    plot(orbit_num, rel_oe_2(:,k));
    grid on;
    xlabel('Orbit Number');
    ylabel(labels{k});
    title(labels{k});
end
