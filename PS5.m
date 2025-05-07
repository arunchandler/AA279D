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
u_TSX_init    = nu_TSX_init + omega_TSX_init;

TSX_init_oe = [a_TSX_init, e_TSX_init, i_TSX_init, ...
               RAAN_TSX_init, omega_TSX_init, M_TSX_init];

% map to TDX via your chosen rel_qns_init
TDX_init_oe = qns2oe(TSX_init_oe, rel_qns_init);

% timing parameters, propagation, etc. (unchanged)
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init^3);
T          = 2*pi/n;
n_orbit    = 15;
tend       = n_orbit*T;
num_points = 1000;
dt         = (tend - tstart)/(num_points-1);
t_grid     = linspace(tstart, tend, num_points).';
t_orbit    = t_grid / T;

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
        % M-C to M-D1
        rel_qns_pre = [0, 0, 0, 300, 0, 400];
        rel_qns_post  = [0, 0, 0, 300, 0, 500];
        scenario_name = 'M-D2 to M-D3';
        
    case 2
        % M-D1 to M-E
        rel_qns_pre = [0, 0, 0, 300, 0, 500];
        rel_qns_post  = [0, 0, 0, 500, 0, 300];
        scenario_name = 'M-D3 to M-D4';
        
    otherwise
        error('Invalid scenario');
end

fprintf('Running scenario "%s": rel_qns_init = [%g %g %g %g %g %g]\n', ...
        'rel_qns_fin = [%g %g %g %g %g %g]\n', ...
        scenario_name, [rel_qns_pre, rel_qns_post]);

TDX_init_oe_pre = qns2oe(TSX_init_oe, rel_qns_pre);

% propagate with your ode4 + J2
[~, TDX_oe_pre] = ode4(@compute_rates_GVE_J2, [tstart, tend/2]', TDX_init_oe_pre,  dt);
TDX_oe_pre(:,6) = wrapTo2Pi(TDX_oe_pre(:,6));

% compute relative OEs and RTN
rel_oe_pre = zeros(num_points/2, 6);
rtn_pre    = zeros(num_points/2, 6);

for idx = 1:num_points/2
    
    a1 = TSX_oe(idx,1); e1 = TSX_oe(idx,2); i1 = TSX_oe(idx,3);
    RAAN1 = TSX_oe(idx,4); omega1 = TSX_oe(idx,5);
    M1 = wrapTo2Pi(TSX_oe(idx,6)); u1 = M1 + omega1;
    
    a2 = TDX_oe_pre(idx,1); e2 = TDX_oe_pre(idx,2); i2 = TDX_oe_pre(idx,3);
    RAAN2 = TDX_oe_pre(idx,4); omega2 = TDX_oe_pre(idx,5);
    M2 = wrapTo2Pi(TDX_oe_pre(idx,6)); u2 = M2 + omega2;
    
    rel_oe_pre(idx,:) = a1 * compute_roes([a1, e1, i1, RAAN1, omega1, M1], ...
                                      [a2, e2, i2, RAAN2, omega2, M2])';
    
    r1 = oe2rv(TSX_oe(idx,:), mu);
    r2 = oe2rv(TDX_oe_pre(idx,:), mu);
    rtn_pre(idx,:) = eci2rtn(r1, r2)';

end

% delta V calculation

Ddiy     = rel_qns_post(6) - rel_qns_pre(6);
Ddey     = rel_qns_post(4) - rel_qns_pre(4);
Ddlambda = rel_qns_post(2) - rel_qns_pre(2);

a_m = TSX_oe(num_points/2,1);

%out of plane portion
uM_op = pi/2;
dvn = Ddiy * n / sin(uM_op)

%in plane portion
dvt = 0.0; %because we dont want sma change, we can set this to 0 and choose dvr to satisfy change in dey
uM_ip1 = 0.0;
uM_ip2 = pi - uM_ip1;
dvr1 = Ddey * n * a_m / (-2 * cos(uM_ip1))
dvr2 = -dvr1

% --- compute total delta‐V ---
DV_op  = abs(dvn);        % out‐of‐plane impulse magnitude
DV_ip1 = abs(dvr1);       % in‐plane impulse at u = 0
DV_ip2 = abs(dvr2);       % in‐plane impulse at u = pi
DV_tot = DV_op + DV_ip1 + DV_ip2;

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

%Post manuever propagation
TDX_init_oe_post = qns2oe(TSX_oe(num_points/2+1,:), rel_qns_post);
[~, TDX_oe_post] = ode4(@compute_rates_GVE_J2, [tstart, tend/2]', TDX_init_oe_post, dt);
TDX_oe_post(:,6) = wrapTo2Pi(TDX_oe_post(:,6));

rel_oe_post = zeros(num_points/2, 6);
rtn_post    = zeros(num_points/2, 6);

for idx = 1:num_points/2
    
    a1 = TSX_oe(num_points/2+idx,1); e1 = TSX_oe(num_points/2+idx,2); i1 = TSX_oe(num_points/2+idx,3);
    RAAN1 = TSX_oe(num_points/2+idx,4); omega1 = TSX_oe(num_points/2+idx,5);
    M1 = wrapTo2Pi(TSX_oe(num_points/2+idx,6)); u1 = M1 + omega1;
    
    a2 = TDX_oe_post(idx,1); e2 = TDX_oe_post(idx,2); i2 = TDX_oe_post(idx,3);
    RAAN2 = TDX_oe_post(idx,4); omega2 = TDX_oe_post(idx,5);
    M2 = wrapTo2Pi(TDX_oe_post(idx,6)); u2 = M2 + omega2;
    
    rel_oe_post(idx,:) = a1 * compute_roes([a1, e1, i1, RAAN1, omega1, M1], ...
                                      [a2, e2, i2, RAAN2, omega2, M2])';
    
    r1 = oe2rv(TSX_oe(num_points/2+idx,:), mu);
    r2 = oe2rv(TDX_oe_post(idx,:), mu);
    rtn_post(idx,:) = eci2rtn(r1, r2)';

end

% RTN frame: pre- and post-maneuver
rR_pre  = rtn_pre(:,1); 
rT_pre  = rtn_pre(:,2); 
rN_pre  = rtn_pre(:,3);
rR_post = rtn_post(:,1);
rT_post = rtn_post(:,2);
rN_post = rtn_post(:,3);

figure;
subplot(1,3,1)
plot(rT_pre,  rR_pre,  'b-', ...
     rT_post, rR_post, 'r--'); grid on; axis equal
xlabel('r_T [m]'); ylabel('r_R [m]');
title('RT plane'); legend('Pre','Post')

subplot(1,3,2)
plot(rN_pre,  rR_pre,  'b-', ...
     rN_post, rR_post, 'r--'); grid on; axis equal
xlabel('r_N [m]'); ylabel('r_R [m]');
title('NR plane'); legend('Pre','Post')

subplot(1,3,3)
plot(rT_pre,  rN_pre,  'b-', ...
     rT_post, rN_post, 'r--'); grid on; axis equal
xlabel('r_T [m]'); ylabel('r_N [m]');
title('TN plane'); legend('Pre','Post')


% Relative orbital elements: pre- and post-maneuver
da_pre      = rel_oe_pre(:,1);
dex_pre     = rel_oe_pre(:,2);
dey_pre     = rel_oe_pre(:,3);
dix_pre     = rel_oe_pre(:,4);
diy_pre     = rel_oe_pre(:,5);
dlambda_pre = rel_oe_pre(:,6);

da_post      = rel_oe_post(:,1);
dex_post     = rel_oe_post(:,2);
dey_post     = rel_oe_post(:,3);
dix_post     = rel_oe_post(:,4);
diy_post     = rel_oe_post(:,5);
dlambda_post = rel_oe_post(:,6);

figure;
subplot(1,3,1)
plot(da_pre,      dlambda_pre,      'b-', ...
     da_post,     dlambda_post,     'r--'); grid on; axis equal
xlabel('a\delta a [m]'); ylabel('a\delta \lambda [m]');
title('a\delta a vs a\delta \lambda'); legend('Pre','Post')

subplot(1,3,2)
plot(dex_pre,     dey_pre,          'b-', ...
     dex_post,    dey_post,         'r--'); grid on; axis equal
xlabel('a\delta e_x [m]'); ylabel('a\delta e_y [m]');
title('a\delta e_x vs a\delta e_y'); legend('Pre','Post')

subplot(1,3,3)
plot(dix_pre,     diy_pre,          'b-', ...
     dix_post,    diy_post,         'r--'); grid on; axis equal
xlabel('a\delta i_x [m]'); ylabel('a\delta i_y [m]');
title('a\delta i_x vs a\delta i_y'); legend('Pre','Post')