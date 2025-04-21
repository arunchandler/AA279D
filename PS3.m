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
i_TDX_init = deg2rad(97.4454);
RAAN_TDX_init = deg2rad(351.0106);
omega_TDX_init = deg2rad(100.5043);
M_TDX_init = deg2rad(12.35936); % added 0.0001

global n

n = sqrt(mu/a_TSX_init^3);
T = 2*pi/n;
t_out_orbit = t_vec / T;

% Part b) init rv

%initial cheif rv
rv_TSX_init = oe2rv([a_TSX_init, e_TSX_init, i_TSX_init, RAAN_TSX_init, omega_TSX_init, M_TSX_init], mu);

%initial deputy rv
rv_TDX_init = oe2rv([a_TDX_init, e_TDX_init, i_TDX_init, RAAN_TDX_init, omega_TDX_init, M_TDX_init], mu);

%relativive rv
[rv_rel_init_RTN, ~] = eci2rtn(rv_TSX_init, rv_TDX_init);

%orbital element differeces
a_diff_init = a_TDX_init - a_TSX_init;
i_diff_init = i_TDX_init - i_TSX_init;
e_diff_init = e_TDX_init - e_TSX_init;
RAAN_diff_init = RAAN_TSX_init - RAAN_TDX_init;
omega_diff_init = omega_TSX_init - omega_TDX_init;
M_diff_init = M_TSX_init - M_TDX_init;

% Part c) compute six HCW integration constants

Ks_init = getHCWconstants(rv_rel_init_RTN, a_TSX_init, tstart);

% Part d) propagate state using HCW equations

[t_out, TDX_RTN] = ode4(@compute_rates_rv_HCW_unperturbed, [tstart, tend]', rv_rel_init_RTN, tint);

figure;
subplot(3,1,1)
plot(t_out_orbit, TDX_RTN(:,1));
xlabel('Orbits');
ylabel('R Position [m]');
grid on;

subplot(3,1,2)
plot(t_out_orbit, TDX_RTN(:,2));
xlabel('Orbits');
ylabel('T Position [m]');
grid on;

subplot(3,1,3)
plot(t_out_orbit, TDX_RTN(:,3));
xlabel('Orbits');
ylabel('N Position [m]');
grid on;

figure;
subplot(2,2,1)
plot(TDX_RTN(:,2), TDX_RTN(:,1))
xlabel('T Position [m]')
ylabel('R Position [m]')
grid on

subplot(2,2,2)
plot(TDX_RTN(:,3), TDX_RTN(:,1))
xlabel('N Position [m]')
ylabel('R Position [m]')
grid on

subplot(2,2,3)
plot(TDX_RTN(:,2), TDX_RTN(:,3))
xlabel('T Position [m]')
ylabel('N Position [m]')
grid on

subplot(2,2,4)
plot3(TDX_RTN(:,1), TDX_RTN(:,2), TDX_RTN(:,3))
xlabel('R Position [m]')
ylabel('T Position [m]')
zlabel('N Position [m]')
grid on
axis equal
view(3)

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
omega_TDX_init_2 = deg2rad(100.5043);
M_TDX_init_2 = deg2rad(12.35936); % added 0.0001


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
RAAN_diff_init_2 = RAAN_TSX_init_2 - RAAN_TDX_init_2;
omega_diff_init_2 = omega_TSX_init_2 - omega_TDX_init_2;
M_diff_init_2 = M_TSX_init_2 - M_TDX_init_2;

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
N       = 1000;
t_grid  = linspace(0, n_orbit*T, N).';    % [s]
t_orbit = t_grid / T;

% preallocate RTN arrays
r_RTN = zeros(N,3);
v_RTN = zeros(N,3);

for k = 1:N
    t = t_grid(k);

    % 1) Mean → Eccentric → True anomaly
    M = M0 + n*t;
    E = mean2ecc(M, e, tol);     % your provided solver
    f = ecc2true(E, e);          % your provided converter

    % 2) Build Phi at (f, f0) with the correct tau
    tau = n*t / eta^3;           % ∫df/k^2 = n·t/η^3
    Phi = buildYAphi(a, e, f, f0, tau);

    % 3) normalized state & back to RTN
    X      = Phi * Ks;           % [x̄; ȳ; z̄; x̄'; ȳ'; z̄']
    r_RTN(k,:) = (X(1:3)' * r0);

    % 4) recover velocity
    r_k     = p / (1 + e*cos(f));  
    fdot    = h / r_k^2;
    v_RTN(k,:) = (X(4:6)' * r0) * fdot;
end

figure;
subplot(3,1,1)
plot(t_orbit, r_RTN(:,1)), xlabel('Orbits'), ylabel('R [m]'), grid on
subplot(3,1,2)
plot(t_orbit, r_RTN(:,2)), xlabel('Orbits'), ylabel('T [m]'), grid on
subplot(3,1,3)
plot(t_orbit, r_RTN(:,3)), xlabel('Orbits'), ylabel('N [m]'), grid on

figure;
subplot(2,2,1)
plot(r_RTN(:,2), r_RTN(:,1)), xlabel('T [m]'), ylabel('R [m]'), grid on
subplot(2,2,2)
plot(r_RTN(:,3), r_RTN(:,1)), xlabel('N [m]'), ylabel('R [m]'), grid on
subplot(2,2,3)
plot(r_RTN(:,2), r_RTN(:,3)), xlabel('T [m]'), ylabel('N [m]'), grid on
subplot(2,2,4)
plot3(r_RTN(:,1), r_RTN(:,2), r_RTN(:,3))
xlabel('R [m]'), ylabel('T [m]'), zlabel('N [m]'),
axis equal, view(3), grid on

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

% time grid
n_orbit = 15;
T       = 2*pi/sqrt(mu/a_TSX_init_2^3);
t_grid  = linspace(0, n_orbit*T, 1000).';

[r_RTN_lin, v_RTN_lin] = propagateLinearEcc(roe_qns, chief_oe, t_grid, tol);

figure('Name','Relative Position Comparison','NumberTitle','off');

% R‐component
subplot(3,1,1)
plot(t_orbit, r_RTN(:,1),'b--', 'DisplayName','YA'); hold on;
plot(t_orbit, r_RTN_lin(:,1),'r-',  'DisplayName','Linear');
ylabel('R [m]'); legend('Location','best'); grid on

% T‐component
subplot(3,1,2)
plot(t_orbit, r_RTN(:,2),'b--'); hold on;
plot(t_orbit, r_RTN_lin(:,2),'r-');
ylabel('T [m]'); grid on

% N‐component
subplot(3,1,3)
plot(t_orbit, r_RTN(:,3),'b--'); hold on;
plot(t_orbit, r_RTN_lin(:,3),'r-');
xlabel('Orbits'); ylabel('N [m]'); grid on


figure('Name','Relative Velocity Comparison','NumberTitle','off');

% Ṙ‐component
subplot(3,1,1)
plot(t_orbit, v_RTN(:,1),'b--','DisplayName','YA'); hold on;
plot(t_orbit, v_RTN_lin(:,1),'r-','DisplayName','Linear');
ylabel('Ṙ [m/s]'); legend('Location','best'); grid on

% Ṫ‐component
subplot(3,1,2)
plot(t_orbit, v_RTN(:,2),'b--'); hold on;
plot(t_orbit, v_RTN_lin(:,2),'r-');
ylabel('Ṫ [m/s]'); grid on

% Ṅ‐component
subplot(3,1,3)
plot(t_orbit, v_RTN(:,3),'b--'); hold on;
plot(t_orbit, v_RTN_lin(:,3),'r-');
xlabel('Orbits'); ylabel('Ṅ [m/s]'); grid on

figure('Name','Relative Paths','NumberTitle','off');

subplot(2,2,1)
plot(r_RTN(:,2), r_RTN(:,1), 'b--'); hold on;
plot(r_RTN_lin(:,2), r_RTN_lin(:,1), 'r-')
xlabel('T [m]'), ylabel('R [m]'), grid on

subplot(2,2,2)
plot(r_RTN(:,3), r_RTN(:,1), 'b--'); hold on; 
plot(r_RTN_lin(:,3), r_RTN_lin(:,1), 'r-')
xlabel('N [m]'), ylabel('R [m]'), grid on

subplot(2,2,3)
plot(r_RTN(:,2), r_RTN(:,3), 'b--'); hold on;
plot(r_RTN_lin(:,2), r_RTN_lin(:,3), 'r-')
xlabel('T [m]'), ylabel('N [m]'), grid on

subplot(2,2,4)
plot3(r_RTN(:,1), r_RTN(:,2), r_RTN(:,3), 'b--'); hold on;
plot3(r_RTN_lin(:,1), r_RTN_lin(:,2), r_RTN_lin(:,3), 'r-')
xlabel('R [m]'), ylabel('T [m]'), zlabel('N [m]')
axis equal, grid on, view(3)
