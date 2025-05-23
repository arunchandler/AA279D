clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

tstart = 0.0;
tint = 1.0;
tend = s_d*1.0;
t_vec = [tstart:tint:tend];
num_points = floor((tend-tstart)/tint) + 1;

%% Problem 1
% Part a) init oe

%initial cheif elements
a_TSX_init = 6886536.686;  % m
e_TSX_init = 0.0001264;
i_TSX_init = deg2rad(97.4453);
RAAN_TSX_init = deg2rad(351.0108);
omega_TSX_init = deg2rad(101.2452);
M_TSX_init = deg2rad(11.6520);
nu_TSX_init = mean2true(M_TSX_init, e_TSX_init,tol);
u_TSX_init = nu_TSX_init + omega_TSX_init;

%initial deputy elements
a_TDX_init =  6886536.686;
e_TDX_init = 0.0001269;
i_TDX_init = deg2rad(97.4554);
RAAN_TDX_init = deg2rad(351.0106);
omega_TDX_init = deg2rad(-100.5043);
M_TDX_init = deg2rad(201.086+12.34926);

global n

n = sqrt(mu/a_TSX_init^3);
T = 2*pi/n;
t_out_orbit = t_vec / T;

% Part b) init rv

%initial cheif rv
rv_TSX_init = oe2rv([a_TSX_init, e_TSX_init, i_TSX_init, RAAN_TSX_init, omega_TSX_init, M_TSX_init], mu)

%initial deputy rv
rv_TDX_init = oe2rv([a_TDX_init, e_TDX_init, i_TDX_init, RAAN_TDX_init, omega_TDX_init, M_TDX_init], mu)

%relativive rv
[rv_rel_init_RTN, ~] = eci2rtn(rv_TSX_init, rv_TDX_init);
disp(rv_rel_init_RTN)

%orbital element differeces
a_diff_init = a_TDX_init - a_TSX_init;
i_diff_init = i_TDX_init - i_TSX_init;
e_diff_init = e_TDX_init - e_TSX_init;
RAAN_diff_init = RAAN_TDX_init - RAAN_TSX_init;
omega_diff_init = omega_TDX_init - omega_TSX_init;
M_diff_init = M_TDX_init - M_TSX_init;

% Display results
fprintf('Semi-major axis difference: %.6f m\n', a_diff_init);
fprintf('Inclination difference: %.6f deg\n', rad2deg(i_diff_init));
fprintf('Eccentricity difference: %.6f\n', e_diff_init);
fprintf('RAAN difference: %.6f deg\n', rad2deg(RAAN_diff_init));
fprintf('Argument of perigee difference: %.6f deg\n', rad2deg(omega_diff_init));
fprintf('Mean anomaly difference: %.6f deg\n', rad2deg(M_diff_init));

% Part c) compute six HCW integration constants

Ks_init = getHCWconstants(rv_rel_init_RTN, a_TSX_init, tstart)

% Part d) propagate state using HCW equations

TDX_RTN_HCW = zeros(num_points, 6);

%[t_out, TDX_RTN] = ode4(@compute_rates_rv_HCW_unperturbed, [tstart, tend]', rv_rel_init_RTN, tint);

for idx = 1:num_points
    Phi = buildHCWphi(a_TSX_init, t_vec(idx));
    TDX_RTN_HCW(idx, :) = Phi * Ks_init;
end

%% Problem 1 d) HCW propagation – positions vs orbits
figure('Name','HCW RTN Position vs Orbits','NumberTitle','off');
subplot(3,1,1)
  plot(t_out_orbit, TDX_RTN_HCW(:,1));
  xlabel('Orbits'); ylabel('R [m]'); grid on;
subplot(3,1,2)
  plot(t_out_orbit, TDX_RTN_HCW(:,2));
  xlabel('Orbits'); ylabel('T [m]'); grid on;
subplot(3,1,3)
  plot(t_out_orbit, TDX_RTN_HCW(:,3));
  xlabel('Orbits'); ylabel('N [m]'); grid on;

%% Problem 1 d) HCW propagation – 2D & 3D RTN trajectories
figure('Name','HCW RTN Trajectories','NumberTitle','off');
subplot(2,2,1)
  plot(TDX_RTN_HCW(:,2), TDX_RTN_HCW(:,1))
  xlabel('T [m]'); ylabel('R [m]'); grid on;
subplot(2,2,2)
  plot(TDX_RTN_HCW(:,3), TDX_RTN_HCW(:,1))
  xlabel('N [m]'); ylabel('R [m]'); grid on;
subplot(2,2,3)
  plot(TDX_RTN_HCW(:,2), TDX_RTN_HCW(:,3))
  xlabel('T [m]'); ylabel('N [m]'); grid on;
subplot(2,2,4)
  plot3(TDX_RTN_HCW(:,1), TDX_RTN_HCW(:,2), TDX_RTN_HCW(:,3))
  xlabel('R [m]'); ylabel('T [m]'); zlabel('N [m]');
  axis equal; grid on; view(3);

%% Problem 2

% a) init oe for TH

%initial cheif elements
a_TSX_init_2 = 6886536.686;  % m
e_TSX_init_2 = 0.1; % changed to 0.1 for cheif orbit
i_TSX_init_2 = deg2rad(97.4453);
RAAN_TSX_init_2 = deg2rad(351.0108);
omega_TSX_init_2 = deg2rad(101.2452);
M_TSX_init_2 = deg2rad(11.6520);
nu_TSX_init_2 = mean2true(M_TSX_init, e_TSX_init,tol);
u_TSX_init_2 = nu_TSX_init + omega_TSX_init;

%initial deputy elements
a_TDX_init_2 =  6886536.686; % still same semi-major axis 
e_TDX_init_2 = 0.1 + 5e-6; % small offset in eccentricity but close enough to keep distance ratio low
i_TDX_init_2 = deg2rad(97.4454);
RAAN_TDX_init_2 = deg2rad(351.0106);
omega_TDX_init_2 = deg2rad(-100.5043);
M_TDX_init_2 = deg2rad(201.086+12.35936); % added 0.0001

%initial cheif rv
rv_TSX_init_2 = oe2rv([a_TSX_init_2, e_TSX_init_2, i_TSX_init_2, RAAN_TSX_init_2, omega_TSX_init_2, M_TSX_init_2], mu);

%initial deputy rv
rv_TDX_init_2 = oe2rv([a_TDX_init_2, e_TDX_init_2, i_TDX_init_2, RAAN_TDX_init_2, omega_TDX_init_2, M_TDX_init_2], mu);

%relativive rv
[rv_rel_init_RTN_2, ~] = eci2rtn(rv_TSX_init_2, rv_TDX_init_2);

%orbital element differeces
a_diff_init_2 = a_TDX_init_2 - a_TSX_init_2;
i_diff_init_2 = i_TDX_init_2 - i_TSX_init_2;
e_diff_init_2 = e_TDX_init_2 - e_TSX_init_2;
RAAN_diff_init_2 = RAAN_TDX_init_2 - RAAN_TSX_init_2;
omega_diff_init_2 = omega_TDX_init_2 - omega_TSX_init_2;
M_diff_init_2 = M_TDX_init_2 - M_TSX_init_2;

r_peri = a_TSX_init_2*(1 - e_TSX_init_2);

ro_init = norm(rv_TSX_init_2(1:3) - rv_TDX_init_2(1:3));

p  = a_TSX_init_2 * (1 - e_TSX_init_2^2);
r0  = p / (1 + e_TSX_init_2 * cos(nu_TSX_init_2));
ratio = ro_init / r_peri; 

fprintf('\nRelevant Initial Orbital Elements for TH Formulation\n');
fprintf('---------------------------------------------------\n');
fprintf('Initial Minimum Distance      : %.15f m\n', r_peri);
fprintf('Initial Relative Distance       : %.15f m\n', ro_init);
fprintf('Initial Minimum Distance Other Way       : %.15f m\n', r0);
fprintf('Initial Relative Distance Other Way       : %.15f m\n', ro_init/r0);
fprintf('Initial Distance Ratio: %.15f \n', ratio);
fprintf('Initial Chief Eccentricity: %.15f \n', e_TSX_init_2);

a    = a_TSX_init_2;
e    = e_TSX_init_2;
f0   = nu_TSX_init_2;
M0   = M_TSX_init_2;
n    = sqrt(mu/a^3);
p    = a*(1 - e^2);
h    = sqrt(mu*p);
eta  = sqrt(1 - e^2);
r0   = p / (1 + e*cos(f0));

%— b) compute your integration constants once ——
Ks = getYAconstants(rv_rel_init_RTN_2, a, e, f0);

fprintf('\nYA Integration Constants:\n');
fprintf('-------------------------\n');
for ii = 1:6
    fprintf('K%d = %12.6e\n', ii, Ks(ii));
end

%— c) set up your time/orbit grid ——
n_orbit = 15;
T       = 2*pi / n;
N       = 10000;
t_grid  = linspace(0, n_orbit*T, N).';    % [s]
t_orbit = t_grid / T;

% preallocate RTN arrays
r_RTN = zeros(N,3);
v_RTN = zeros(N,3);

n_orbit = 15;
T       = 2*pi/n;
N       = 10000;
t_grid  = linspace(0, n_orbit*T, N).';

r_RTN = zeros(N,3);
v_RTN = zeros(N,3);

for k = 1:N
    t = t_grid(k);

    % 1) Mean → Eccentric → True anomaly
    M = M0 + n*t;
    f = mean2true(M, e, tol);

    % 2) Build YA transition (absolute units)
    tau = n*t / eta^3;
    Phi = buildYAphi(a, e, f, f0, tau);

    % 3) Apply to your K‐vector
    X = Phi * Ks;      % now X = [x; y; z; vx; vy; vz] in [m,m,m,m/s,m/s,m/s]

    % 4) Store directly
    r_RTN(k,:) = X(1:3)';    % [m]
    v_RTN(k,:) = X(4:6)';    % [m/s]
end


% d) explain trends (need to write up) 

% e) quasi-nonsingular orbital elements

da = (a_TDX_init_2 - a_TSX_init_2) / a_TSX_init_2;
dlambda = (M_TDX_init_2 + omega_TDX_init_2) - (M_TSX_init_2 + omega_TSX_init_2) + (RAAN_TDX_init_2 - RAAN_TSX_init_2)*cos(i_TSX_init_2);
de_x = e_TDX_init_2*cos(omega_TDX_init_2) - e_TSX_init_2*cos(omega_TSX_init_2);
de_y = e_TDX_init_2*sin(omega_TDX_init_2) - e_TSX_init_2*sin(omega_TSX_init_2);
di_x = i_TDX_init_2 - i_TSX_init_2;
di_y = (RAAN_TDX_init_2 - RAAN_TSX_init_2)*sin(i_TSX_init_2); 


% Print ROEs for June 21, 2010
fprintf('\nRelative Quasi-Nonsingular Orbit Elements:\n');
fprintf('---------------------------------------------------\n');
fprintf('Relative Semi-Major Axis      : %.15f \n', da);
fprintf('Relative Lambda       : %.15f \n', dlambda);
fprintf('Relative Eccentricity Vector X: %.15f \n', de_x);
fprintf('Relative Eccentricity Vector Y: %.15f \n', de_y);
fprintf('Relative Inclination Vector X : %.15f \n', di_x);
fprintf('Relative Inclination Vector Y : %.15f \n', di_y);

roe_qns  = [da; dlambda; de_x; de_y; di_x; di_y];
chief_oe = [a_TSX_init_2; e_TSX_init_2; i_TSX_init_2; ...
            RAAN_TSX_init_2; omega_TSX_init_2; M_TSX_init_2];

% f) geometric linear mapping propagation
% time grid
n_orbit = 15;
T       = 2*pi/sqrt(mu/a_TSX_init_2^3);
t_int   = n_orbit*T/(N-1);
t_start = 0;
t_end = n_orbit*T;
t_grid  = linspace(t_start, t_end, N).';

[r_RTN_lin, v_RTN_lin] = propagateLinearEcc(roe_qns, chief_oe, t_grid, tol);

% h) nonlinear equations propagation

[t_out, TDX_RTN_2] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [t_start, t_end]', [rv_rel_init_RTN_2; rv_TSX_init_2], t_int);

% check rho ratio throughout entire diff eq simulation: 

r_rel = TDX_RTN_2(:,1:3);            % [m] (N×3)

% preallocate
Npts = numel(t_out);
r_chief = zeros(Npts,1);

for k = 1:Npts
    Mk = M0 + n * t_out(k);
    Ek = mean2ecc(Mk, e, tol);
    fk = ecc2true(Ek, e);
    r_chief(k) = p / (1 + e*cos(fk));    % [m]
end

% 3) compute the ratio vector and its average
rho    = vecnorm(r_rel, 2, 2);          % ||relative|| [m]
ratio  = rho ./ r_chief;                % dimensionless
ratio_avg = mean(ratio);

fprintf('\nNonlinear Trajectory Separation Check\n');
fprintf('-------------------------------------\n');
fprintf('Mean ||ρ||/r_c  = %.3e (should be ≲1e–3)\n', ratio_avg);
fprintf('Max  ||ρ||/r_c  = %.3e\n', max(ratio));
fprintf('Min  ||ρ||/r_c  = %.3e\n', min(ratio));



% convert to km / km-s
rYA   = r_RTN       / 1e3;         % [km]
rLIN  = r_RTN_lin   / 1e3;         % [km]
rTR   = TDX_RTN_2(:,1:3) / 1e3;    % [km]

%vTR   = TDX_RTN_2(:,4:6);    % [m/s]


% 1) YA only — Position Planes
figure('Name','YA Only (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rYA(:,2), rYA(:,1), 'r', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('R [km]'); grid on;
subplot(1,3,2)
  plot(rYA(:,3), rYA(:,1), 'r', 'LineWidth', 1.0);
  xlabel('N [km]'); ylabel('R [km]'); grid on;
subplot(1,3,3)
  plot(rYA(:,2), rYA(:,3), 'r', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 1b) YA only — Position 3D
figure('Name','YA Only (3D Pos)','NumberTitle','off');
plot3(rYA(:,1), rYA(:,2), rYA(:,3), 'r', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
axis equal; grid on; view(3);


% 2) YA + Linear — Position Planes
figure('Name','YA vs Linear (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rYA(:,2), rYA(:,1),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,2),rLIN(:,1),'b', 'LineWidth', 1.0); 
  xlabel('T [km]'); ylabel('R [km]'); grid on;
  legend('YA','Geo. Map.','Location','best');
subplot(1,3,2)
  plot(rYA(:,3), rYA(:,1),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,3),rLIN(:,1),'b', 'LineWidth', 1.0); 
  xlabel('N [km]'); ylabel('R [km]'); grid on;
subplot(1,3,3)
  plot(rYA(:,2), rYA(:,3),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0); 
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 2b) YA + Linear — Position 3D
figure('Name','YA vs Linear (3D Pos)','NumberTitle','off');
plot3(rYA(:,1),rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0); hold on;
plot3(rLIN(:,1),rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
legend('YA','Geo. Map.','Location','best'); axis equal; grid on; view(3);


% 3) YA + Linear + Truth — Position Planes
figure('Name','All Three (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rTR(:,2),rTR(:,1),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,2),rYA(:,1),'r', 'LineWidth', 1.0);
  plot(rLIN(:,2),rLIN(:,1),'b', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('R [km]'); grid on;
  legend('Non. Linear Diff. Eq. ','YA','Geo. Map.','Location','best');
subplot(1,3,2)
  plot(rTR(:,3),rTR(:,1),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,3),rYA(:,1),'r', 'LineWidth', 1.0);
  plot(rLIN(:,3),rLIN(:,1),'b', 'LineWidth', 1.0);
  xlabel('N [km]'); ylabel('R [km]');grid on;
subplot(1,3,3)
  plot(rTR(:,2),rTR(:,3),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0);
  plot(rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 3b) YA + Linear + Truth — Position 3D
figure('Name','All Three (3D Pos)','NumberTitle','off');
plot3(rTR(:,1),rTR(:,2),rTR(:,3),'c', 'LineWidth', 1.0); hold on;
plot3(rYA(:,1),rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0);
plot3(rLIN(:,1),rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
legend('Non. Linear Diff. Eq. ','YA','Geo. Map.','Location','best');
axis equal; grid on; view(3);

% i) Problem 2 but with da

%initial cheif elements
a_TSX_init_2 = 6886536.686;  % m
e_TSX_init_2 = 0.0001264;
i_TSX_init_2 = deg2rad(97.4453);
RAAN_TSX_init_2 = deg2rad(351.0108);
omega_TSX_init_2 = deg2rad(101.2452);
M_TSX_init_2 = deg2rad(11.6520);
nu_TSX_init_2 = mean2true(M_TSX_init, e_TSX_init,tol);
u_TSX_init_2 = nu_TSX_init + omega_TSX_init;

%initial deputy elements
a_TDX_init_2 =  6886436.686; % -100
e_TDX_init_2 = 0.0001269;
i_TDX_init_2 = deg2rad(97.4454);
RAAN_TDX_init_2 = deg2rad(351.0106);
omega_TDX_init_2 = deg2rad(-100.5043);
M_TDX_init_2 = deg2rad(201.086+12.35936); % added 0.0001


%initial cheif rv
rv_TSX_init_2 = oe2rv([a_TSX_init_2, e_TSX_init_2, i_TSX_init_2, RAAN_TSX_init_2, omega_TSX_init_2, M_TSX_init_2], mu);

%initial deputy rv
rv_TDX_init_2 = oe2rv([a_TDX_init_2, e_TDX_init_2, i_TDX_init_2, RAAN_TDX_init_2, omega_TDX_init_2, M_TDX_init_2], mu);

%relativive rv
[rv_rel_init_RTN_2, ~] = eci2rtn(rv_TSX_init_2, rv_TDX_init_2);

%orbital element differeces
a_diff_init_2 = a_TDX_init_2 - a_TSX_init_2;
i_diff_init_2 = i_TDX_init_2 - i_TSX_init_2;
e_diff_init_2 = e_TDX_init_2 - e_TSX_init_2;
RAAN_diff_init_2 = RAAN_TDX_init_2 - RAAN_TSX_init_2;
omega_diff_init_2 = omega_TDX_init_2 - omega_TSX_init_2;
M_diff_init_2 = M_TDX_init_2 - M_TSX_init_2;

r_peri = a_TSX_init_2*(1 - e_TSX_init_2);

ro_init = norm(rv_TSX_init_2(1:3) - rv_TDX_init_2(1:3));

p  = a_TSX_init_2 * (1 - e_TSX_init_2^2);
r0  = p / (1 + e_TSX_init_2 * cos(nu_TSX_init_2));
ratio = ro_init / r_peri; 

fprintf('\nPart i-i: Relevant Initial Orbital Elements for TH Formulation \n');
fprintf('---------------------------------------------------\n');
fprintf('Initial Minimum Distance      : %.15f m\n', r_peri);
fprintf('Initial Relative Distance       : %.15f m\n', ro_init);
fprintf('Initial Minimum Distance Other Way       : %.15f m\n', r0);
fprintf('Initial Relative Distance Other Way       : %.15f m\n', ro_init/r0);
fprintf('Initial Distance Ratio: %.15f \n', ratio);
fprintf('Initial Chief Eccentricity: %.15f \n', e_TSX_init_2);

a    = a_TSX_init_2;
e    = e_TSX_init_2;
f0   = nu_TSX_init_2;
M0   = M_TSX_init_2;
n    = sqrt(mu/a^3);
p    = a*(1 - e^2);
h    = sqrt(mu*p);
eta  = sqrt(1 - e^2);
r0   = p / (1 + e*cos(f0));


Ks = getYAconstants(rv_rel_init_RTN_2, a, e, f0);

fprintf('\nPart i-i YA Integration Constants:\n');
fprintf('-------------------------\n');
for ii = 1:6
    fprintf('K%d = %12.6e\n', ii, Ks(ii));
end

n_orbit = 15;
T       = 2*pi/n;
N       = 10000;
t_grid  = linspace(0, n_orbit*T, N).';

r_RTN = zeros(N,3);
v_RTN = zeros(N,3);

for k = 1:N
    t = t_grid(k);

    % 1) Mean → Eccentric → True anomaly
    M = M0 + n*t;
    f = mean2true(M, e, tol);

    % 2) Build YA transition (absolute units)
    tau = n*t / eta^3;
    Phi = buildYAphi(a, e, f, f0, tau);

    % 3) Apply to your K‐vector
    X = Phi * Ks;      % now X = [x; y; z; vx; vy; vz] in [m,m,m,m/s,m/s,m/s]

    % 4) Store directly
    r_RTN(k,:) = X(1:3)';    % [m]
    v_RTN(k,:) = X(4:6)';    % [m/s]
end


da = (a_TDX_init_2 - a_TSX_init_2) / a_TSX_init_2;
dlambda = (M_TDX_init_2 + omega_TDX_init_2) - (M_TSX_init_2 + omega_TSX_init_2) + (RAAN_TDX_init_2 - RAAN_TSX_init_2)*cos(i_TSX_init_2);
de_x = e_TDX_init_2*cos(omega_TDX_init_2) - e_TSX_init_2*cos(omega_TSX_init_2);
de_y = e_TDX_init_2*sin(omega_TDX_init_2) - e_TSX_init_2*sin(omega_TSX_init_2);
di_x = i_TDX_init_2 - i_TSX_init_2;
di_y = (RAAN_TDX_init_2 - RAAN_TSX_init_2)*sin(i_TSX_init_2); 


% Print ROEs for June 21, 2010
fprintf('\nPart i-i Relative Quasi-Nonsingular Orbit Elements:\n');
fprintf('---------------------------------------------------\n');
fprintf('Relative Semi-Major Axis      : %.15f \n', da);
fprintf('Relative Lambda       : %.15f \n', dlambda);
fprintf('Relative Eccentricity Vector X: %.15f \n', de_x);
fprintf('Relative Eccentricity Vector Y: %.15f \n', de_y);
fprintf('Relative Inclination Vector X : %.15f \n', di_x);
fprintf('Relative Inclination Vector Y : %.15f \n', di_y);

roe_qns  = [da; dlambda; de_x; de_y; di_x; di_y];
chief_oe = [a_TSX_init_2; e_TSX_init_2; i_TSX_init_2; ...
            RAAN_TSX_init_2; omega_TSX_init_2; M_TSX_init_2];

% time grid
n_orbit = 15;
T       = 2*pi/sqrt(mu/a_TSX_init_2^3);
t_int   = n_orbit*T/(N-1);
t_start = 0;
t_end = n_orbit*T;
t_grid  = linspace(t_start, t_end, N).';

[r_RTN_lin, v_RTN_lin] = propagateLinearEcc(roe_qns, chief_oe, t_grid, tol);



[t_out, TDX_RTN_2] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [t_start, t_end]', [rv_rel_init_RTN_2; rv_TSX_init_2], t_int);

% h) nonlinear equations propagation

[t_out, TDX_RTN_2] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [t_start, t_end]', [rv_rel_init_RTN_2; rv_TSX_init_2], t_int);

% check rho ratio throughout entire diff eq simulation: 

r_rel = TDX_RTN_2(:,1:3);            % [m] (N×3)

% preallocate
Npts = numel(t_out);
r_chief = zeros(Npts,1);

for k = 1:Npts
    Mk = M0 + n * t_out(k);
    Ek = mean2ecc(Mk, e, tol);
    fk = ecc2true(Ek, e);
    r_chief(k) = p / (1 + e*cos(fk));    % [m]
end

% 3) compute the ratio vector and its average
rho    = vecnorm(r_rel, 2, 2);          % ||relative|| [m]
ratio  = rho ./ r_chief;                % dimensionless
ratio_avg = mean(ratio);

fprintf('\nPart i-i Nonlinear Trajectory Separation Check\n');
fprintf('-------------------------------------\n');
fprintf('Mean ||ρ||/r_c  = %.3e (should be ≲1e–3)\n', ratio_avg);
fprintf('Max  ||ρ||/r_c  = %.3e\n', max(ratio));
fprintf('Min  ||ρ||/r_c  = %.3e\n', min(ratio));

 
% 1) YA only — Position Planes
figure('Name','δa = –100 m:YA Only (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rYA(:,2), rYA(:,1), 'r', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('R [km]'); grid on;
subplot(1,3,2)
  plot(rYA(:,3), rYA(:,1), 'r', 'LineWidth', 1.0);
  xlabel('N [km]'); ylabel('R [km]'); grid on;
subplot(1,3,3)
  plot(rYA(:,2), rYA(:,3), 'r', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 1b) YA only — Position 3D
figure('Name','δa = –100 m:YA Only (3D Pos)','NumberTitle','off');
plot3(rYA(:,1), rYA(:,2), rYA(:,3), 'r', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
axis equal; grid on; view(3);



% 2) YA + Linear — Position Planes
figure('Name','δa = –100 m: YA vs Linear (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rYA(:,2), rYA(:,1),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,2),rLIN(:,1),'b', 'LineWidth', 1.0); 
  xlabel('T [km]'); ylabel('R [km]'); grid on;
  legend('YA','Geo. Map.','Location','best');
subplot(1,3,2)
  plot(rYA(:,3), rYA(:,1),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,3),rLIN(:,1),'b', 'LineWidth', 1.0); 
  xlabel('N [km]'); ylabel('R [km]'); grid on;
subplot(1,3,3)
  plot(rYA(:,2), rYA(:,3),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0); 
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 2b) YA + Linear — Position 3D
figure('Name','δa = –100 m: YA vs Linear (3D Pos)','NumberTitle','off');
plot3(rYA(:,1),rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0); hold on;
plot3(rLIN(:,1),rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
legend('YA','Geo. Map.','Location','best'); axis equal; grid on; view(3);




% 3) YA + Linear + Truth — Position Planes
figure('Name','δa = –100 m: All Three (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rTR(:,2),rTR(:,1),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,2),rYA(:,1),'r', 'LineWidth', 1.0);
  plot(rLIN(:,2),rLIN(:,1),'b', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('R [km]'); grid on;
  legend('Non. Linear Diff. Eq. ','YA','Geo. Map.','Location','best');
subplot(1,3,2)
  plot(rTR(:,3),rTR(:,1),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,3),rYA(:,1),'r', 'LineWidth', 1.0);
  plot(rLIN(:,3),rLIN(:,1),'b', 'LineWidth', 1.0);
  xlabel('N [km]'); ylabel('R [km]');grid on;
subplot(1,3,3)
  plot(rTR(:,2),rTR(:,3),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0);
  plot(rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 3b) YA + Linear + Truth — Position 3D
figure('Name','δa = –100 m: All Three (3D Pos)','NumberTitle','off');
plot3(rTR(:,1),rTR(:,2),rTR(:,3),'c', 'LineWidth', 1.0); hold on;
plot3(rYA(:,1),rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0);
plot3(rLIN(:,1),rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
legend('Non. Linear Diff. Eq. ','YA','Geo. Map.','Location','best');
axis equal; grid on; view(3);

% ii) Problem 2 but with e > 0.5

%initial cheif elements
a_TSX_init_2 = 6886536.686;  % m
e_TSX_init_2 = 0.6; % changed to 0.1 for cheif orbit
i_TSX_init_2 = deg2rad(97.4453);
RAAN_TSX_init_2 = deg2rad(351.0108);
omega_TSX_init_2 = deg2rad(101.2452);
M_TSX_init_2 = deg2rad(11.6520);
nu_TSX_init_2 = mean2true(M_TSX_init, e_TSX_init,tol);
u_TSX_init_2 = nu_TSX_init + omega_TSX_init;

%initial deputy elements
a_TDX_init_2 =  6886536.686; % still same semi-major axis 
e_TDX_init_2 = 0.6;
i_TDX_init_2 = deg2rad(97.4454);
RAAN_TDX_init_2 = deg2rad(351.0106);
omega_TDX_init_2 = deg2rad(-100.5043);
M_TDX_init_2 = deg2rad(201.086+12.35936); % added 0.0001


%initial cheif rv
rv_TSX_init_2 = oe2rv([a_TSX_init_2, e_TSX_init_2, i_TSX_init_2, RAAN_TSX_init_2, omega_TSX_init_2, M_TSX_init_2], mu);

%initial deputy rv
rv_TDX_init_2 = oe2rv([a_TDX_init_2, e_TDX_init_2, i_TDX_init_2, RAAN_TDX_init_2, omega_TDX_init_2, M_TDX_init_2], mu);

%relativive rv
[rv_rel_init_RTN_2, ~] = eci2rtn(rv_TSX_init_2, rv_TDX_init_2);

%orbital element differeces
a_diff_init_2 = a_TDX_init_2 - a_TSX_init_2;
i_diff_init_2 = i_TDX_init_2 - i_TSX_init_2;
e_diff_init_2 = e_TDX_init_2 - e_TSX_init_2;
RAAN_diff_init_2 = RAAN_TDX_init_2 - RAAN_TSX_init_2;
omega_diff_init_2 = omega_TDX_init_2 - omega_TSX_init_2;
M_diff_init_2 = M_TDX_init_2 - M_TSX_init_2;

r_peri = a_TSX_init_2*(1 - e_TSX_init_2);

ro_init = norm(rv_TSX_init_2(1:3) - rv_TDX_init_2(1:3));

p  = a_TSX_init_2 * (1 - e_TSX_init_2^2);
r0  = p / (1 + e_TSX_init_2 * cos(nu_TSX_init_2));
ratio = ro_init / r_peri; 

fprintf('\nPart i-iiRelevant Initial Orbital Elements for TH Formulation\n');
fprintf('---------------------------------------------------\n');
fprintf('Initial Minimum Distance      : %.15f m\n', r_peri);
fprintf('Initial Relative Distance       : %.15f m\n', ro_init);
fprintf('Initial Minimum Distance Other Way       : %.15f m\n', r0);
fprintf('Initial Relative Distance Other Way       : %.15f m\n', ro_init/r0);
fprintf('Initial Distance Ratio: %.15f \n', ratio);
fprintf('Initial Chief Eccentricity: %.15f \n', e_TSX_init_2);

a    = a_TSX_init_2;
e    = e_TSX_init_2;
f0   = nu_TSX_init_2;
M0   = M_TSX_init_2;
n    = sqrt(mu/a^3);
p    = a*(1 - e^2);
h    = sqrt(mu*p);
eta  = sqrt(1 - e^2);
r0   = p / (1 + e*cos(f0));

Ks = getYAconstants(rv_rel_init_RTN_2, a, e, f0);

fprintf('\nPart i-ii YA Integration Constants:\n');
fprintf('-------------------------\n');
for ii = 1:6
    fprintf('K%d = %12.6e\n', ii, Ks(ii));
end

n_orbit = 15;
T       = 2*pi/n;
N       = 10000;
t_grid  = linspace(0, n_orbit*T, N).';

r_RTN = zeros(N,3);
v_RTN = zeros(N,3);

for k = 1:N
    t = t_grid(k);

    % 1) Mean → Eccentric → True anomaly
    M = M0 + n*t;
    f = mean2true(M, e, tol);

    % 2) Build YA transition (absolute units)
    tau = n*t / eta^3;
    Phi = buildYAphi(a, e, f, f0, tau);

    % 3) Apply to your K‐vector
    X = Phi * Ks;      % now X = [x; y; z; vx; vy; vz] in [m,m,m,m/s,m/s,m/s]

    % 4) Store directly
    r_RTN(k,:) = X(1:3)';    % [m]
    v_RTN(k,:) = X(4:6)';    % [m/s]
end

da = (a_TDX_init_2 - a_TSX_init_2) / a_TSX_init_2;
dlambda = (M_TDX_init_2 + omega_TDX_init_2) - (M_TSX_init_2 + omega_TSX_init_2) + (RAAN_TDX_init_2 - RAAN_TSX_init_2)*cos(i_TSX_init_2);
de_x = e_TDX_init_2*cos(omega_TDX_init_2) - e_TSX_init_2*cos(omega_TSX_init_2);
de_y = e_TDX_init_2*sin(omega_TDX_init_2) - e_TSX_init_2*sin(omega_TSX_init_2);
di_x = i_TDX_init_2 - i_TSX_init_2;
di_y = (RAAN_TDX_init_2 - RAAN_TSX_init_2)*sin(i_TSX_init_2); 


% Print ROEs for June 21, 2010
fprintf('\n Part i-ii nRelative Quasi-Nonsingular Orbit Elements:\n');
fprintf('---------------------------------------------------\n');
fprintf('Relative Semi-Major Axis      : %.15f \n', da);
fprintf('Relative Lambda       : %.15f \n', dlambda);
fprintf('Relative Eccentricity Vector X: %.15f \n', de_x);
fprintf('Relative Eccentricity Vector Y: %.15f \n', de_y);
fprintf('Relative Inclination Vector X : %.15f \n', di_x);
fprintf('Relative Inclination Vector Y : %.15f \n', di_y);

roe_qns  = [da; dlambda; de_x; de_y; di_x; di_y];
chief_oe = [a_TSX_init_2; e_TSX_init_2; i_TSX_init_2; ...
            RAAN_TSX_init_2; omega_TSX_init_2; M_TSX_init_2];

% time grid
n_orbit = 15;
T       = 2*pi/sqrt(mu/a_TSX_init_2^3);
t_int   = n_orbit*T/(N-1);
t_start = 0;
t_end = n_orbit*T;
t_grid  = linspace(t_start, t_end, N).';

[r_RTN_lin, v_RTN_lin] = propagateLinearEcc(roe_qns, chief_oe, t_grid, tol);


% h) nonlinear equations propagation

[t_out, TDX_RTN_2] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [t_start, t_end]', [rv_rel_init_RTN_2; rv_TSX_init_2], t_int);

% check rho ratio throughout entire diff eq simulation: 

r_rel = TDX_RTN_2(:,1:3);            % [m] (N×3)

% preallocate
Npts = numel(t_out);
r_chief = zeros(Npts,1);

for k = 1:Npts
    Mk = M0 + n * t_out(k);
    Ek = mean2ecc(Mk, e, tol);
    fk = ecc2true(Ek, e);
    r_chief(k) = p / (1 + e*cos(fk));    % [m]
end

% 3) compute the ratio vector and its average
rho    = vecnorm(r_rel, 2, 2);          % ||relative|| [m]
ratio  = rho ./ r_chief;                % dimensionless
ratio_avg = mean(ratio);

fprintf('\nNonlinear Trajectory Separation Check\n');
fprintf('-------------------------------------\n');
fprintf('Mean ||ρ||/r_c  = %.3e (should be ≲1e–3)\n', ratio_avg);
fprintf('Max  ||ρ||/r_c  = %.3e\n', max(ratio));
fprintf('Min  ||ρ||/r_c  = %.3e\n', min(ratio));

 
% convert to km / km-s
rYA   = r_RTN       / 1e3;         % [km]
rLIN  = r_RTN_lin   / 1e3;         % [km]
rTR   = TDX_RTN_2(:,1:3) / 1e3;    % [km]

%vTR   = TDX_RTN_2(:,4:6);    % [m/s]


% 1) YA only — Position Planes
figure('Name','e = 0.6: YA Only (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rYA(:,2), rYA(:,1), 'r', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('R [km]'); grid on;
subplot(1,3,2)
  plot(rYA(:,3), rYA(:,1), 'r', 'LineWidth', 1.0);
  xlabel('N [km]'); ylabel('R [km]'); grid on;
subplot(1,3,3)
  plot(rYA(:,2), rYA(:,3), 'r', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 1b) YA only — Position 3D
figure('Name','e = 0.6: YA Only (3D Pos)','NumberTitle','off');
plot3(rYA(:,1), rYA(:,2), rYA(:,3), 'r', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
axis equal; grid on; view(3);



% 2) YA + Linear — Position Planes
figure('Name','e = 0.6: YA vs Linear (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rYA(:,2), rYA(:,1),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,2),rLIN(:,1),'b', 'LineWidth', 1.0); 
  xlabel('T [km]'); ylabel('R [km]'); grid on;
  legend('YA','Geo. Map.','Location','best');
subplot(1,3,2)
  plot(rYA(:,3), rYA(:,1),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,3),rLIN(:,1),'b', 'LineWidth', 1.0); 
  xlabel('N [km]'); ylabel('R [km]'); grid on;
subplot(1,3,3)
  plot(rYA(:,2), rYA(:,3),'r', 'LineWidth', 1.0); hold on;
  plot(rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0); 
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 2b) YA + Linear — Position 3D
figure('Name','e = 0.6: YA vs Linear (3D Pos)','NumberTitle','off');
plot3(rYA(:,1),rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0); hold on;
plot3(rLIN(:,1),rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
legend('YA','Geo. Map.','Location','best'); axis equal; grid on; view(3);




% 3) YA + Linear + Truth — Position Planes
figure('Name','e = 0.6: All Three (Pos)','NumberTitle','off');
subplot(1,3,1)
  plot(rTR(:,2),rTR(:,1),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,2),rYA(:,1),'r', 'LineWidth', 1.0);
  plot(rLIN(:,2),rLIN(:,1),'b', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('R [km]'); grid on;
  legend('Non. Linear Diff. Eq. ','YA','Geo. Map.','Location','best');
subplot(1,3,2)
  plot(rTR(:,3),rTR(:,1),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,3),rYA(:,1),'r', 'LineWidth', 1.0);
  plot(rLIN(:,3),rLIN(:,1),'b', 'LineWidth', 1.0);
  xlabel('N [km]'); ylabel('R [km]');grid on;
subplot(1,3,3)
  plot(rTR(:,2),rTR(:,3),'c', 'LineWidth', 1.0); hold on;
  plot(rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0);
  plot(rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
  xlabel('T [km]'); ylabel('N [km]'); grid on;

% 3b) YA + Linear + Truth — Position 3D
figure('Name','e = 0.6: All Three (3D Pos)','NumberTitle','off');
plot3(rTR(:,1),rTR(:,2),rTR(:,3),'c', 'LineWidth', 1.0); hold on;
plot3(rYA(:,1),rYA(:,2),rYA(:,3),'r', 'LineWidth', 1.0);
plot3(rLIN(:,1),rLIN(:,2),rLIN(:,3),'b', 'LineWidth', 1.0);
xlabel('R [km]'); ylabel('T [km]'); zlabel('N [km]');
legend('Non. Linear Diff. Eq. ','YA','Geo. Map.','Location','best');
axis equal; grid on; view(3);