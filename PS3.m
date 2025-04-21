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
