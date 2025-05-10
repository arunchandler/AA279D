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


%% === CONTROL SECTION: Initialize ===
roe_target = [0, 0, 0, 300, 0, -1000];


%% === 1) DETECTION & Δv CALCULATION ===
a_chief = a_TSX_init;                     % chief SMA  [m]
n       = sqrt(mu/a_chief^3);             % mean motion [rad s-1]

% ---- nominal windows (in metres) --------------------------------------
delta_e_max_m = 2.0;                      % half-window ‖a·δe‖
delta_i_max_m = 2.0;                      % half-window ‖a·δi‖

% ---- convert the “design” numbers & ROE targets to dimension-less -----
e_nom_hat = roe_target(3:4)/a_chief;      % δe_nom   (-)
i_nom_hat = roe_target(5:6)/a_chief;      % δi_nom   (-)
de_max    = delta_e_max_m/a_chief;        % ‖δe‖ window (-)
di_max    = delta_i_max_m/a_chief;        % ‖δi‖ window (-)

angle_thresh_deg = 7;                     % trigger geometry
T_half = pi/n;                            % ½-orbit time [s]
v_chief = sqrt(mu / a_chief);
% helper: 2-D rotation about +z
Rz = @(phi)[cos(phi) -sin(phi); sin(phi) cos(phi)];

maneuver_log = [];                        % save first cycle only

for idx = 1:num_points
    
    % ---------- current relative elements (dimension-less) -------------
    roe_hat  = rel_oe(idx,:)'/a_chief;    % [δa δλ δex δey δix δiy]'
    da_hat   = roe_hat(1);                % δa
    de_hat   = roe_hat(3:4);              % δe  (x,y)
    di_hat   = roe_hat(5:6);              % δi  (x,y)
    
    % ---------- chief argument of latitude (for timing) ----------------
    M_c  = wrapTo2Pi(TSX_oe(idx,6));
    nu_c = mean2true(M_c, TSX_oe(idx,2), tol);
    u_c  = wrapTo2Pi(nu_c + TSX_oe(idx,5));
    
    % ===================================================================
    % 1)   Post-manoeuvre targets  (Eqs. 2.58–2.62)  — *all dimension-less*
    % ===================================================================
    gamma   = 0.75*J2*(Re/a_chief)^2;
    phi_dot = 1.5*gamma*(5*cos(i_TSX_init)^2 - 1);   % sign only is needed
    dphi    = sign(phi_dot)*asin(de_max/ norm(e_nom_hat));
    
    de_man_hat = (Rz(dphi)*e_nom_hat.').';                  % Eq. (2.58)
    di_man_hat = [ i_nom_hat(1);
                   i_nom_hat(2) - sign(i_nom_hat(1))*di_max ]; % Eq. (2.61)
    
    % ===================================================================
    % 2)   Corrections still required & pulse longitudes
    % ===================================================================
    de_req = de_man_hat - de_hat;          % δe that must be fixed
    di_req = di_man_hat - di_hat;          % δi that must be fixed
    
    u1 = wrapTo2Pi(atan2(de_req(2), de_req(1)));   % first tangential
    u2 = mod(u1 + pi, 2*pi);                       % second tangential
    un = wrapTo2Pi(atan2(di_req(2), di_req(1)));   % cross-track
    
    %% 3) Choose the maintenance-cycle length  Δt  [s]
    %    (time from this 3-burn set to the *next* one)
    Delta_t = 2*T_half;          % ← 1 full orbit; tune to your mission
    
    %% 3-a)  secular J2 + drag contributions to δu  (Eq. 2.69)
    du_J2 = -12*gamma*sin(2*i_TSX_init)*di_hat(1)*n*Delta_t;   % J2 term
    du_D  = 0;                                                 % drag term
                                                    %  (set to 0 if neglected)
    
    %% 3-b)  build the bracket of Eq. (2.71)
    delta_u      = roe_hat(2);          % current along-track error  ̂δu
    delta_u_nom  = 0;                   % centre of the window (your target)
    
    bracket =  3*de_max                     ...   % 3 δe_max
            +  da_hat                       ...   % + δa
            - (4/(3*pi)) * ( delta_u - delta_u_nom ...
                             + du_J2 + du_D );     % – 4/3π(…).
    
    %% 3-c)  δa^{man} from Eq. (2.71)
    den     = 2*n*Delta_t - pi;          % denominator  2 n Δt − π
    da_man  = -pi/den * bracket;         % ***dimension-less***
    
    %% 4)  Δv magnitudes (Eqs. 2.63–2.64) with the *dimension-less* variables
    dvt1 = (n*a_chief/4) * ( (da_man - da_hat) + norm(de_req) );
    dvt2 = (n*a_chief/4) * ( (da_man - da_hat) - norm(de_req) );
    dvn  =  n*a_chief     *   norm(di_req);
    
    % ===================================================================
    % 5)   Trigger once per cycle (unchanged)
    % ===================================================================
    sep = rad2deg(acos(dot(de_hat, di_hat)/(norm(de_hat)*norm(di_hat))));
    sep = min(sep, 180 - sep);
    
    if sep > angle_thresh_deg && isempty(maneuver_log)
        
        maneuver_log = struct( ...
            'idx',   idx,        't_det', t_grid(idx), ...
            'dvt1',  dvt1,       'dvt2',  dvt2,        ...
            'dvn',   dvn,        'u1',    u1,          ...
            'u2',    u2,         'un',    un,          ...
            'u_c',   u_c                         );
        
        fprintf('[%6.2f orbits] trigger at idx=%d\n', t_grid(idx)/T, idx);
        fprintf('  Δv_T1 = %.4f m/s @ u1 = %.2f rad\n', dvt1, u1);
        fprintf('  Δv_T2 = %.4f m/s @ u2 = %.2f rad\n', dvt2, u2);
        fprintf('  Δv_N  = %.4f m/s @ un = %.2f rad (sep = %.2f°)\n\n', ...
                dvn, un, sep);
    end
end

%% === 2) APPLY ALL THREE MANEUVERS ===
if ~isempty(maneuver_log)
    det = maneuver_log;

    % 1) First along‑track burn at u1
    t1       = det.t_det + wrapTo2Pi(det.u1 - det.u_c)/n;
    [~, i1]  = min(abs(t_grid - t1));
    oe1_m    = TDX_oe(i1,:);  nu1 = mean2true(oe1_m(6),oe1_m(2),tol);
    oe1_t    = [oe1_m(1:5), nu1];
    rv1      = oe2rv(oe1_t,mu);
    R1       = rv1(1:3)/norm(rv1(1:3));
    N1       = cross(rv1(1:3),rv1(4:6)); N1=N1/norm(N1);
    T1       = cross(N1,R1);
    dv1_vec  = det.dvt1 * T1;
    v1_new   = rv1(4:6) + dv1_vec;
    oe1_f    = rv2oe([rv1(1:3);v1_new], mu);
    M1_f     = true2mean(oe1_f(6),oe1_f(2));
    oe1_new  = [oe1_f(1:5), M1_f];

    % 2) Half‑orbit propagate to second T‐burn
    t2     = t1 + T_half;
    [~, i2s] = min(abs(t_grid - t1));
    [~, i2e] = min(abs(t_grid - t2));
    dt_fine  = (t2 - t1)/(i2e-i2s+1);
    [~, oe_mid] = ode4(@compute_rates_GVE_J2, [t1 t2]', oe1_new', dt_fine);

    % 3) Second along‑track burn at u2
    oe2_m   = oe_mid(end,:); nu2 = mean2true(oe2_m(6),oe2_m(2),tol);
    oe2_t   = [oe2_m(1:5), nu2];
    rv2     = oe2rv(oe2_t, mu);
    R2      = rv2(1:3)/norm(rv2(1:3));
    N2      = cross(rv2(1:3),rv2(4:6)); N2=N2/norm(N2);
    T2      = cross(N2,R2);
    dv2_vec = det.dvt2 * T2;
    v2_new  = rv2(4:6) + dv2_vec;
    oe2_f   = rv2oe([rv2(1:3);v2_new], mu);
    M2_f    = true2mean(oe2_f(6),oe2_f(2));
    oe2_new = [oe2_f(1:5), M2_f];

    % 4) Cross‑track burn at un
    t3     = det.t_det + wrapTo2Pi(det.un - det.u_c)/n;
    [~, i3] = min(abs(t_grid - t3));
    oe3_m   = TDX_oe(i3,:); nu3 = mean2true(oe3_m(6),oe3_m(2),tol);
    oe3_t   = [oe3_m(1:5), nu3];
    rv3     = oe2rv(oe3_t,mu);
    N3      = cross(rv3(1:3),rv3(4:6)); N3=N3/norm(N3);
    v3_new  = rv3(4:6) + det.dvn * N3;
    oe3_f   = rv2oe([rv3(1:3);v3_new], mu);
    M3_f    = true2mean(oe3_f(6),oe3_f(2));
    oe3_new = [oe3_f(1:5), M3_f];

    % 5) Final propagate
    [~, oe_end] = ode4(@compute_rates_GVE_J2, [t3 tend]', oe3_new', dt);

    % 6) Stitch TDX back
    TDX_oe = [ TDX_oe(1:i1-1,:)
               oe1_new
               oe_mid(2:end,:)
               oe2_new
               TDX_oe(i2e+1:i3-1,:)
               oe3_new
               oe_end(2:end,:) ];

    % 7) Print final ROE & e/i sep:
    rel_fin = TSX_oe(i3,1)*compute_roes(TSX_oe(i3,:),TDX_oe(i3,:))';
    de  = norm(rel_fin(3:4)); di = norm(rel_fin(5:6));
    sep = rad2deg(acos(dot(rel_fin(3:4),rel_fin(5:6))/(de*di)));
    fprintf('\n=== Post‑burn @ idx=%d ===\n', i3);
    fprintf(' Δv_T1=%.4f, Δv_T2=%.4f, Δv_N=%.4f m/s\n',det.dvt1,det.dvt2,det.dvn);
    fprintf(' ROE=[%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f]\n', rel_fin);
    fprintf(' ||δe||=%6.2f m, ||δi||=%6.2f m, sep=%6.2f°\n\n',de,di,sep);
end