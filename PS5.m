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

