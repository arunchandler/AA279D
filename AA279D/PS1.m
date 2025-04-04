clc; clear; close all;

format long g;

tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;

%% Problem 1
fprintf('\n ---------PROBLEM 1------- \n');
% calculating semi-major axis
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400;

mean_motions = [15.19150740 15.19152819 15.20699603 15.19152132];

n_rad = mean_motions .* (2 * pi) ./ s_d;

a_m = (mu ./ (n_rad .^ 2)).^(1/3);

fprintf('Mean Motion (rev/day)  |  Semi-Major Axis (m)\n');
fprintf('----------------------------------------------\n');
for i = 1:length(mean_motions)
    fprintf('%10.6f         |  %f m\n', mean_motions(i), a_m(i));
end

% calculating relative orbital elements
% Epoch 1: June 21, 2010
a_TSX_1 = 6886542.969;   % m
i_TSX_1 = deg2rad(97.4435); % rad
e_TSX_1 = 0.0001655;
RAAN_TSX_1 = deg2rad(179.2734); % rad
omega_TSX_1 = deg2rad(84.9587);
M_TSX_1 = deg2rad(342.8671);

a_TDX_1 = 6881866.120;  % m
i_TDX_1 = deg2rad(97.4512); % rad
e_TDX_1 = 0.0011015;
RAAN_TDX_1 = deg2rad(179.2734); % rad
omega_TDX_1 = deg2rad(285.9794);
M_TDX_1 = deg2rad(74.0206);

% Epoch 2: Dec 12, 2010
a_TSX_2 = 6886536.686;  % m
i_TSX_2 = deg2rad(97.4453); % rad
e_TSX_2 = 0.0001264;
RAAN_TSX_2 = deg2rad(351.0108); % rad
omega_TSX_2 = deg2rad(101.2452);
M_TSX_2 = deg2rad(11.6520);

a_TDX_2 = 6886538.762;  % m
i_TDX_2 = deg2rad(97.4454); % rad
e_TDX_2 = 0.0001269;
RAAN_TDX_2 = deg2rad(351.0106); % rad
omega_TDX_2 = deg2rad(100.5043);
M_TDX_2 = deg2rad(12.3926);

nu_TSX_1 = mean2true(M_TSX_1, e_TSX_1,tol);
nu_TDX_1 = mean2true(M_TDX_1, e_TDX_1,tol);
u_TSX_1 = nu_TSX_1 + omega_TSX_1;
u_TDX_1 = nu_TDX_1 + omega_TDX_1;

nu_TSX_2 = mean2true(M_TSX_2, e_TSX_2,tol);
nu_TDX_2 = mean2true(M_TDX_2, e_TDX_2,tol);
u_TSX_2 = nu_TSX_2 + omega_TSX_2;
u_TDX_2 = nu_TDX_2 + omega_TDX_2;

roe_1 = compute_roes(a_TSX_1, i_TSX_1, e_TSX_1, RAAN_TSX_1, u_TSX_1, a_TDX_1, i_TDX_1, e_TDX_1, RAAN_TDX_1, u_TDX_1);
roe_2 = compute_roes(a_TSX_2, i_TSX_2, e_TSX_2, RAAN_TSX_2, u_TSX_2, a_TDX_2, i_TDX_2, e_TDX_2, RAAN_TDX_2, u_TDX_2);

% Print ROEs for June 21, 2010
fprintf('\nRelative Orbital Elements (ROEs) for June 21, 2010:\n');
fprintf('---------------------------------------------------\n');
fprintf('Relative Semi-Major Axis      : %.15f m\n', roe_1(1));
fprintf('Relative Mean Longitude       : %.15f rad\n', roe_1(2));
fprintf('Relative Eccentricity Vector X: %.15f \n', roe_1(3));
fprintf('Relative Eccentricity Vector Y: %.15f \n', roe_1(4));
fprintf('Relative Inclination Vector X : %.15f \n', roe_1(5));
fprintf('Relative Inclination Vector Y : %.15f \n', roe_1(6));

% Print ROEs for Dec 12, 2010
fprintf('\nRelative Orbital Elements (ROEs) for Dec 12, 2010:\n');
fprintf('--------------------------------------------------\n');
fprintf('Relative Semi-Major Axis      : %.15f m\n', roe_2(1));
fprintf('Relative Mean Longitude       : %.15f rad\n', roe_2(2));
fprintf('Relative Eccentricity Vector X: %.15f \n', roe_2(3));
fprintf('Relative Eccentricity Vector Y: %.15f \n', roe_2(4));
fprintf('Relative Inclination Vector X : %.15f \n', roe_2(5));
fprintf('Relative Inclination Vector Y : %.15f \n', roe_2(6));


%% Problem 2
fprintf('\n -------PROBLEM 2------- \n');
% part a) Define initial conditions - already done above. Using orbital
% elements at epoch 1

% part b) Treat initial conditions as osculating quantities and compute
% corresponding initial position and velocity in ECI
init_koe = [a_TDX_1, e_TDX_1,i_TDX_1, RAAN_TDX_1, omega_TDX_1, M_TDX_1];
osc_elements = mean2osc(init_koe);

fprintf('Semi-major axis (a): %.6f m\n', osc_elements(1));
fprintf('Eccentricity (e): %.6f\n', osc_elements(2));
fprintf('Inclination (I): %.6f deg\n', rad2deg(osc_elements(3)));
fprintf('RAAN (Ω): %.6f deg\n', rad2deg(osc_elements(4)));
fprintf('Argument of Perigee (ω): %.6f deg\n', rad2deg(osc_elements(5)));
fprintf('Mean Anomaly (M): %.6f deg\n', rad2deg(osc_elements(6)));

init_TDX1_rv_ECI = oe2rv(osc_elements, mu);
r_eci = init_TDX1_rv_ECI(1:3); % Position (m)
v_eci = init_TDX1_rv_ECI(4:6); % Velocity (m/s)

fprintf('\nInitial TDX1 State Vector in ECI Coordinates:\n');
fprintf('----------------------------------------------\n');
fprintf('Position (ECI) [m]:\n');
fprintf('  x: %.15f\n', r_eci(1));
fprintf('  y: %.15f\n', r_eci(2));
fprintf('  z: %.15f\n', r_eci(3));

fprintf('\nVelocity (ECI) [m/s]:\n');
fprintf('  vx: %.15f\n', v_eci(1));
fprintf('  vy: %.15f\n', v_eci(2));
fprintf('  vz: %.15f\n', v_eci(3));

%part c) - numerically propagate state (rv) including and excluidng J2
tstart = 0.0;
tint = 10.0;
tend = s_d * 5;

%unperturbed
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
%[t_out, TDX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TDX1_rv_ECI, tint);
[t_out, TDX_rv_out_unperturbed] = ode45(@compute_rates_rv_unperturbed, [tstart:tint:tend]', init_TDX1_rv_ECI, options);

figure;
plot3(TDX_rv_out_unperturbed(:,1), TDX_rv_out_unperturbed(:,2), TDX_rv_out_unperturbed(:,3));
hold on;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('TDX Unperturbed Trajectory in ECI Frame');
grid on;
axis equal;
[X, Y, Z] = sphere(50);  % Generate sphere points
surf(Re * X, Re * Y, Re * Z, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

%perturbed
%[t_out, TDX_rv_out_perturbed] = ode4(@compute_rates_rv_perturbed, [tstart, tend]', init_TDX1_rv_ECI, tint);
[t_out, TDX_rv_out_perturbed] = ode45(@compute_rates_rv_perturbed, [tstart:tint:tend]', init_TDX1_rv_ECI, options);

figure;
plot3(TDX_rv_out_perturbed(:,1), TDX_rv_out_perturbed(:,2), TDX_rv_out_perturbed(:,3));
hold on;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
title('TDX Perturbed Trajectory in ECI Frame');
grid on;
axis equal;
surf(Re * X, Re * Y, Re * Z, 'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

%part d) analytical Kepplerian propagation
a = osc_elements(1);
M0 = osc_elements(6);
n = sqrt(mu/a^3);
M_arr = wrapTo2Pi(M0 + n*t_out);

n_steps = length(t_out);
rv_eci_matrix = zeros(n_steps, 6);

for i = 1:n_steps
    osc_elements(6) = M_arr(i);
    rv_eci = oe2rv(osc_elements, mu);
    rv_eci_matrix(i, :) = rv_eci';
end

position_differences = rv_eci_matrix(:, 1:3) - TDX_rv_out_unperturbed(:, 1:3);
velocity_differences = rv_eci_matrix(:, 4:6) - TDX_rv_out_unperturbed(:, 4:6);

figure;
% Position differences
subplot(2,1,1);
plot(t_out, vecnorm(position_differences, 2, 2), 'b');
xlabel('Time (s)');
ylabel('Position Difference (m)');
title('Position Difference Over Time');
grid on;
% Velocity differences
subplot(2,1,2);
plot(t_out, vecnorm(velocity_differences, 2, 2), 'r');
xlabel('Time (s)');
ylabel('Velocity Difference (m/s)');
title('Velocity Difference Over Time');
grid on;

%part e) compute koe, ecc vector, angular momentum vector, specific mechanical energy from rv for all TDX out
num_states = size(TDX_rv_out_unperturbed, 1);
params = zeros(num_states, 13);

for k = 1:num_states
    params(k, :) = rv2oe(TDX_rv_out_unperturbed(k, :), mu);
end

a = params(:, 1);
e = params(:, 2);
i = params(:, 3);
RAAN = params(:, 4);
omega = params(:, 5);
nu = params(:, 6);
energy = params(:, 7);
h = params(:, 8:10);
e_vec = params(:, 11:13);

% 1. Plot Semi-major axis (a)
figure;
plot(1:num_states, a, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Semi-major axis (a) [km]');
title('Semi-major axis vs Time');

% 2. Plot Eccentricity (e)
figure;
plot(1:num_states, e, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Eccentricity (e)');
title('Eccentricity vs Time');

% 3. Plot Inclination (i)
figure;
plot(1:num_states, i, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Inclination (i) [deg]');
title('Inclination vs Time');

% 4. Plot RAAN (Right Ascension of Ascending Node)
figure;
plot(1:num_states, RAAN, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('RAAN [deg]');
title('RAAN vs Time');

% 5. Plot Argument of Periapsis (omega)
figure;
plot(1:num_states, omega, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Argument of Periapsis (ω) [deg]');
title('Argument of Periapsis vs Time');

% 6. Plot True Anomaly (nu)
figure;
plot(1:num_states, nu, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('True Anomaly (ν) [deg]');
title('True Anomaly vs Time');

%Energy
figure;
plot(1:num_states, energy, 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Energy (J)');
title('Energy vs Time');

%Angular Momentum (h)
figure;
subplot(3, 1, 1);
plot(1:num_states, h(:, 1), 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('h_x');
title('Angular Momentum (h_x) vs Time');

subplot(3, 1, 2);
plot(1:num_states, h(:, 2), 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('h_y');
title('Angular Momentum (h_y) vs Time');

subplot(3, 1, 3);
plot(1:num_states, h(:, 3), 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('h_z');
title('Angular Momentum (h_z) vs Time');

%Eccentricity Vector (e_vec)
figure;
subplot(3, 1, 1);
plot(1:num_states, e_vec(:, 1), 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('e_x');
title('Eccentricity Vector (e_x) vs Time');

subplot(3, 1, 2);
plot(1:num_states, e_vec(:, 2), 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('e_y');
title('Eccentricity Vector (e_y) vs Time');