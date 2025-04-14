clc; clear; close all;

addpath('mean_osc');

format long g;

global tol Re J2 mu s_d
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
tint = 1.0;
tend = s_d*0.5;

[init_TDX_rv_RTN,~] = eci2rtn(init_TSX_rv_ECI, init_TDX_rv_ECI);
r0_init = init_TSX_rv_ECI;

init_state = [init_TDX_rv_RTN; r0_init];

[t_out, TDX_rv_out_unperturbed_nonlin] = ode4(@compute_rates_rv_rel_unperturbed, [tstart, tend]', init_state, tint);

figure;
sgtitle('Relative Position & Velocity');
subplot(2,3,1)
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,1));
xlabel('Time [s]');
ylabel('R Position [m]');
grid on;

subplot(2,3,2)
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,2));
xlabel('Time [s]');
ylabel('T Position [m]');
grid on;

subplot(2,3,3)
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,3));
xlabel('Time [s]');
ylabel('N Position [m]');
grid on;

subplot(2, 3, 4)
plot(t_out, TDX_rv_out_unperturbed_nonlin(:, 4));
xlabel('Time [s]');
ylabel('R Velocity [m/s]');
grid on;

subplot(2, 3, 5)
plot(t_out, TDX_rv_out_unperturbed_nonlin(:, 5));
xlabel('Time [s]');
ylabel('T Velocity [m/s]');
grid on;

subplot(2, 3, 6)
plot(t_out, TDX_rv_out_unperturbed_nonlin(:, 6));
xlabel('Time [s]');
ylabel('N Velocity [m/s]');
grid on;

figure;
plot3(TDX_rv_out_unperturbed_nonlin(:,1), TDX_rv_out_unperturbed_nonlin(:,2), TDX_rv_out_unperturbed_nonlin(:,3));
hold on;
scatter3(0, 0, 0, 100, 'Marker', '.');
grid on;
axis equal;
xlabel('R Position [m]');
ylabel('T Position [m]');
zlabel('N Position [m]');
legend('TDX Position', 'TSX Position')
title('3D Relative Orbit of TDX in TSX''s RTN Frame');


%% Part c) compute relative orbit with fundamental diff eq

[t_out, TSX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TSX_rv_ECI, tint);
[t_out, TDX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TDX_rv_ECI, tint);

num_points = size(TSX_rv_out_unperturbed, 1);
TDX_rv_RTN = zeros(num_points, 6);

for idx = 1:num_points
    TSX_rv = TSX_rv_out_unperturbed(idx, :)';
    TDX_rv = TDX_rv_out_unperturbed(idx, :)';
    [TDX_rv_RTN(idx, :), ~] = eci2rtn(TSX_rv, TDX_rv);
end

figure;
sgtitle('Relative Position & Velocity Comparison');
subplot(2,3,1)
plot(t_out, TDX_rv_RTN(:,1));
hold on;
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,1));
xlabel('Time [s]');
ylabel('R Position [m]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,2)
plot(t_out, TDX_rv_RTN(:,2));
hold on;
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,2));
xlabel('Time [s]');
ylabel('T Position [m]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,3)
plot(t_out, TDX_rv_RTN(:,3));
hold on;
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,3));
xlabel('Time [s]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
ylabel('N Position [m]');
grid on;

subplot(2,3,4)
plot(t_out, TDX_rv_RTN(:,4));
hold on;
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,4));
xlabel('Time [s]');
ylabel('R Velocity [m/s]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,5)
plot(t_out, TDX_rv_RTN(:,5));
hold on;
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,5));
xlabel('Time [s]');
ylabel('T Velocity [m/s]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,6)
plot(t_out, TDX_rv_RTN(:,6));
hold on;
plot(t_out, TDX_rv_out_unperturbed_nonlin(:,6));
xlabel('Time [s]');
ylabel('N Velocity [m/s]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

figure;
plot3(TDX_rv_RTN(:,1), TDX_rv_RTN(:,2), TDX_rv_RTN(:,3));
hold on;
%plot3(TDX_rv_out_unperturbed_nonlin(:,1), TDX_rv_out_unperturbed_nonlin(:,2), TDX_rv_out_unperturbed_nonlin(:,3));
scatter3(0, 0, 0, 100, 'Marker', '.');
grid on;
xlabel('R Position [m]');
ylabel('T Position [m]');
zlabel('N Position [m]');
axis equal;
legend('TDX Position (fundamental diff eq.)', 'TDX Position (relative nonlin. eq.)', 'TSX Position')
title('3D Relative Orbit of TDX in TSX''s RTN Frame');

%% Part d) compare b and c

err_RTN = TDX_rv_out_unperturbed_nonlin(:, 1:6) - TDX_rv_RTN;

figure;
subplot(2,3,1)
plot(t_out, err_RTN(:,1));
xlabel('Time [s]');
ylabel('Error in R Position [m]');
grid on;

subplot(2,3,2)
plot(t_out, err_RTN(:,2));
xlabel('Time [s]');
ylabel('Error in T Position [m]');
grid on;

subplot(2,3,3)
plot(t_out, err_RTN(:,3));
xlabel('Time [s]');
ylabel('Error in N Position [m]');
grid on;

subplot(2,3,4)
plot(t_out, err_RTN(:,4));
xlabel('Time [s]');
ylabel('Error in R Velocity [m/s]');
grid on;

subplot(2,3,5)
plot(t_out, err_RTN(:,5));
xlabel('Time [s]');
ylabel('Error in T Velocity [m/s]');
grid on;

subplot(2,3,6)
plot(t_out, err_RTN(:,6));
xlabel('Time [s]');
ylabel('Error in N Velocity [m/s]');
grid on;

sgtitle('Error Comparison: Unperturbed Nonlinear Eq. vs. Differential ECI');


%% Part e) most fuel-efficient impulse maneuver by deputy - try integrating analytically the GVEs over a dv

%% Part f) Apply maneuver by repeating b and adding discontinuity in the inertial velocity of deputy