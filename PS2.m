clc; clear; close all;

addpath('mean_osc');

format long g;

tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

%% Part a) 

a_TSX_2 = 6886536.686;  % m
i_TSX_2 = deg2rad(97.4453);
e_TSX_2 = 0.0001264;
RAAN_TSX_2 = deg2rad(351.0108);
omega_TSX_2 = deg2rad(101.2452);
M_TSX_2 = deg2rad(11.6520);
nu_TSX_2 = mean2true(M_TSX_2, e_TSX_2,tol);
u_TSX_2 = nu_TSX_2 + omega_TSX_2;

a_TDX_2 = 6886536.686;  % changed so its exactly the same as stated on the pset
i_TDX_2 = deg2rad(97.4454);
e_TDX_2 = 0.0001269;
RAAN_TDX_2 = deg2rad(351.0106);
omega_TDX_2 = deg2rad(100.5043);
M_TDX_2 = deg2rad(12.3926);

init_TSX_koe = [a_TSX_2, e_TSX_2, i_TSX_2, RAAN_TSX_2, omega_TSX_2, M_TSX_2];
init_TDX_koe = [a_TDX_2, e_TDX_2, i_TDX_2, RAAN_TDX_2, omega_TDX_2, M_TDX_2];

init_TSX_rv_ECI = oe2rv(init_TSX_koe, mu);
init_TDX_rv_ECI = oe2rv(init_TDX_koe, mu);

%% Part b) numerical integration of nonlinear equations of relative motion
tstart = 0.0;
tint = 50.0;
tend = s_d*1;

init_rho = init_TDX_rv_ECI - init_TSX_rv_ECI;


%% Part c) compute relative orbit with fundamental diff eq

[t_out, TSX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TSX_rv_ECI, tint);
[t_out, TDX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TDX_rv_ECI, tint);

diff_rv_ECI = (TSX_rv_out_unperturbed - TDX_rv_out_unperturbed)'; % 6 x n

num_points = size(TSX_rv_out_unperturbed, 1);
TDX_r_RTN = zeros(num_points, 3);

for idx = 1:num_points
    TSX_rv = TSX_rv_out_unperturbed(idx, :)';
    [~, R_rtn2eci_TSX] = eci2rtn(TSX_rv);
    R_eci2rtn_TSX = R_rtn2eci_TSX';
    TDX_r_RTN(idx, :) = R_eci2rtn_TSX * diff_rv_ECI(1:3, idx);
end

figure;
subplot(3,1,1)
plot(t_out, TDX_r_RTN(:,1));
xlabel('Time [s]');
ylabel('R Position [m]');
title('Relative Orbit of TDX in TSX''s RTN Frame');
grid on;

subplot(3,1,2)
plot(t_out, TDX_r_RTN(:,2));
xlabel('Time [s]');
ylabel('T Position [m]');
grid on;

subplot(3,1,3)
plot(t_out, TDX_r_RTN(:,3));
xlabel('Time [s]');
ylabel('N Position [m]');
grid on;

figure;
plot3(TDX_r_RTN(:,1), TDX_r_RTN(:,2), TDX_r_RTN(:,3), 'LineWidth', 1.5);
grid on;
xlabel('R Position [m]');
ylabel('T Position [m]');
zlabel('N Position [m]');
title('3D Relative Orbit of TDX in TSX''s RTN Frame');

%% Part d) compare b and c

%% Part e) most fuel-efficient impulse maneuver by deputy - try integrating analytically the GVEs over a dv

%% Part f) Apply maneuver by repeating b and adding discontinuity in the inertial velocity of deputy