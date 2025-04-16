clc; clear; close all;

addpath('mean_osc');

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

T = 2 * pi ./ n_rad;  % Orbital period in seconds for each mean motion
fprintf('\nOrbital Periods (seconds):\n');
for i = 1:length(T)
    fprintf('Period %d: %.2f s\n', i, T(i));
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
%osc_elements = mean2osc(init_koe);
osc_elements = init_koe;

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
tint = 1.0;
tend = s_d*1;
num_points = floor((tend-tstart)/tint) + 1;
t_vec = [tstart:tint:tend];

t_end_extended = s_d*11;
t_vec_extended = [tstart:tint:t_end_extended];


n = sqrt(mu/a_TSX_1^3);
T = 2*pi/n;
fprintf('Orbit Period %.f\n', T);
t_out_orbit = t_vec / T;

t_out_orbit_ext = t_vec_extended / T;

%unperturbed
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, TDX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TDX1_rv_ECI, tint);
%[t_out, TDX_rv_out_unperturbed] = ode113(@compute_rates_rv_unperturbed, [tstart:tint:tend]', init_TDX1_rv_ECI, options);

figure;
tiledlayout(3,2);
titles = {'X [m]', 'Y [m]', 'Z [m]', 'V_x [m/s]', 'V_y [m/s]', 'V_z [m/s]'};
for i = 1:6
    nexttile;
    plot(t_out_orbit, TDX_rv_out_unperturbed(:,i), 'b');
    xlabel('Orbits');
    ylabel(titles{i});
    grid on;
end

%perturbed
[t_out_ext, TDX_rv_out_perturbed] = ode4(@compute_rates_rv_perturbed, [tstart, t_end_extended]', init_TDX1_rv_ECI, tint);
%[t_out, TDX_rv_out_perturbed] = ode113(@compute_rates_rv_perturbed, [tstart:tint:tend]', init_TDX1_rv_ECI, options);

figure;
tiledlayout(3,2);
for i = 1:6
    nexttile;
    plot(t_out_orbit_ext, TDX_rv_out_perturbed(:,i), 'r');
    xlabel('Orbits');
    ylabel(titles{i});
    grid on;
end




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

pos_error_rtn = zeros(n_steps, 3);
vel_error_rtn = zeros(n_steps, 3);

for i = 1:n_steps
    r_ref = TDX_rv_out_unperturbed(i, 1:3)';
    v_ref = TDX_rv_out_unperturbed(i, 4:6)';

    rv_true_eci = rv_eci_matrix(i, :)';

    [delta_rv_rtn, R_eci2rtn] = eci2rtn(rv_true_eci, [r_ref;v_ref]);

    pos_error_rtn(i, :) = delta_rv_rtn(1:3)';
    vel_error_rtn(i, :) = delta_rv_rtn(4:6)';
end

figure;
hold on;
grid on;
axis equal;

% Numerical (unperturbed) trajectory
plot3(TDX_rv_out_unperturbed(:,1), TDX_rv_out_unperturbed(:,2), TDX_rv_out_unperturbed(:,3), ...
    'b');

% Analytical Keplerian propagation
plot3(rv_eci_matrix(:,1), rv_eci_matrix(:,2), rv_eci_matrix(:,3), ...
    'r--');

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');

legend('Numerical Propagation (Unperturbed)', 'Analytical Keplerian');

figure;

% R Position Error
subplot(2,3,1);
plot(t_out_orbit, pos_error_rtn(:,1), 'r');
xlabel('Orbits');
ylabel('R Position Error [m]');
grid on;

% T Position Error
subplot(2,3,2);
plot(t_out_orbit, pos_error_rtn(:,2), 'g');
xlabel('Orbits');
ylabel('T Position Error [m]');
grid on;

% N Position Error
subplot(2,3,3);
plot(t_out_orbit, pos_error_rtn(:,3), 'b');
xlabel('Orbits');
ylabel('N Position Error [m]');
grid on;

% vR Velocity Error
subplot(2,3,4);
plot(t_out_orbit, vel_error_rtn(:,1), 'r');
xlabel('Orbits');
ylabel('R Velocity Error [m/s]');
grid on;

% vT Velocity Error
subplot(2,3,5);
plot(t_out_orbit, vel_error_rtn(:,2), 'g');
xlabel('Orbits');
ylabel('T Velocity Error [m/s]');
grid on;

% vN Velocity Error
subplot(2,3,6);
plot(t_out_orbit, vel_error_rtn(:,3), 'b');
xlabel('Orbits');
ylabel('N Velocity Error [m/s]');
grid on;



%part e) compute koe, ecc vector, angular momentum vector, specific mechanical energy from rv for all TDX out
num_states = size(TDX_rv_out_unperturbed, 1);
params_perturbed = zeros(num_states, 12);
params_unperturbed = zeros(num_states, 12);

for k = 1:num_states
    params_perturbed(k, :) = rv2oe(TDX_rv_out_perturbed(k, :), mu);
    params_unperturbed(k, :) = rv2oe(TDX_rv_out_unperturbed(k, :), mu);
end

% Extract parameters
a_p = params_perturbed(:, 1);
e_p = params_perturbed(:, 2);
i_p = params_perturbed(:, 3);
RAAN_p = params_perturbed(:, 4);
omega_p = params_perturbed(:, 5);
nu_p = params_perturbed(:, 6);
energy_p = params_perturbed(:, 7);
h_p = params_perturbed(:, 8:10);
e_vec_p = params_perturbed(:, 11:12);

a_u = params_unperturbed(:, 1);
e_u = params_unperturbed(:, 2);
i_u = params_unperturbed(:, 3);
RAAN_u = params_unperturbed(:, 4);
omega_u = params_unperturbed(:, 5);
nu_u = params_unperturbed(:, 6);
energy_u = params_unperturbed(:, 7);
h_u = params_unperturbed(:, 8:10);
e_vec_u = params_unperturbed(:, 11:12);

time_vec = [tstart:tint:tend];

% 1. Semi-major axis
figure;
plot(t_out_orbit, a_p, 'b', t_out_orbit, a_u, 'b--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('Semi-major axis (a) [m]');
%title('Semi-major Axis vs Time');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([6.85e6 6.9e6]);


% 2. Eccentricity
figure;
plot(t_out_orbit, e_p, 'r', t_out_orbit, e_u, 'r--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('Eccentricity');
%title('Eccentricity vs Time');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([0 0.003]);


% 3. Inclination
figure;
plot(t_out_orbit, i_p, 'g', t_out_orbit, i_u, 'g--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('Inclination [rad]');
%title('Inclination vs Time');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([1.7 1.702]);

% 4. RAAN
figure;
plot(t_out_orbit, RAAN_p, 'm', t_out_orbit, RAAN_u, 'm--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('RAAN [rad]');
%title('RAAN vs Time');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([3 3.2]);

% 5. Argument of Periapsis
figure;
plot(t_out_orbit, omega_p, 'c', t_out_orbit, omega_u, 'c--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('Argument of Periapsis [rad]');
%title('Argument of Periapsis vs Time t<10000');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([0 2*pi]);

% 6. True Anomaly
figure;
plot(t_out_orbit, nu_p, 'k', t_out_orbit, nu_u, 'k--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('True Anomaly [rad]');
%title('True Anomaly vs Time');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([0 2*pi]);

% 7. Specific Mechanical Energy 
figure;
plot(t_out_orbit, energy_p, 'b', t_out_orbit, energy_u, 'b--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('Specific Mechanical Energy [J/kg]');
%title('Specific Mechanical Energy vs Time');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([-2.91e7 -2.89e7]);

% 8. Angular Momentum Vector
figure;
subplot(3,1,1);
plot(t_out_orbit, h_p(:,1), 'r', t_out_orbit, h_u(:,1), 'r--', 'LineWidth', 1.5);
ylabel('h_x'); %title('Angular Momentum Vector Components');
legend('Perturbed', 'Unperturbed'); grid on;

subplot(3,1,2);
plot(t_out_orbit, h_p(:,2), 'g', t_out_orbit, h_u(:,2), 'g--', 'LineWidth', 1.5);
ylabel('h_y'); legend('Perturbed', 'Unperturbed'); grid on;

subplot(3,1,3);
plot(t_out_orbit, h_p(:,3), 'b', t_out_orbit, h_u(:,3), 'b--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('h_z'); legend('Perturbed', 'Unperturbed'); grid on;
ylim([-7e9 -6e9]);

% 9. Eccentricity Vector Components
figure;
subplot(2,1,1);
plot(t_out_orbit, e_vec_p(:,1), 'r', t_out_orbit, e_vec_u(:,1), 'r--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('e_x'); %title('Eccentricity Vector Components');
legend('Perturbed', 'Unperturbed'); grid on;
ylim([-0.0005 0.0015]);

subplot(2,1,2);
plot(t_out_orbit, e_vec_p(:,2), 'g', t_out_orbit, e_vec_u(:,2), 'g--', 'LineWidth', 1.5);
xlabel('Orbits'); ylabel('e_y'); legend('Perturbed', 'Unperturbed'); grid on;
ylim([-0.0015 0.0015]);



% part f) propagating orbital elements from J2 effects

u_p = nu_p + omega_p;

e_x_p = e_p.*cos(omega_p);
e_y_p = e_p.*sin(omega_p);

TDX_osc_oe = [a_p, e_x_p, e_y_p, i_p, RAAN_p, u_p];

init_osc_elements = zeros(1,6);
init_osc_elements(1) = init_koe(1); %semi-major axis
init_osc_elements(2) = init_koe(2)*cos(init_koe(5)); %ex
init_osc_elements(3) = init_koe(2)*sin(init_koe(5)); %ey
init_osc_elements(4:5) = init_koe(3:4); %i, RAAN
init_osc_elements(6) = mean2true(init_koe(6), init_koe(2), tol) + init_koe(5); %u

[t_out, TDX_osc_oe_propagated] = ode4(@compute_rates_osc_elements, [tstart, tend]', init_osc_elements, tint);

ex_prop = TDX_osc_oe_propagated(:,2);
ey_prop = TDX_osc_oe_propagated(:,3);
RAAN_prop = TDX_osc_oe_propagated(:,5);
u_prop = TDX_osc_oe_propagated(:,6);

%short period
figure;

subplot(4,1,1)
plot(t_out_orbit, e_x_p, 'b', t_out_orbit, ex_prop, 'r--')
ylabel('e_x')
legend('State Propagation', 'OE Propagation')
grid on

subplot(4,1,2)
plot(t_out_orbit, e_y_p, 'b', t_out_orbit, ey_prop, 'r--')
ylabel('e_y')
legend('State Propagation', 'OE Propagation')
grid on

subplot(4,1,3)
plot(t_out_orbit, RAAN_p, 'b', t_out_orbit, RAAN_prop, 'r--')
ylabel('RAAN [rad]')
legend('State Propagation', 'OE Propagation')
grid on

subplot(4,1,4)
plot(t_out_orbit, u_p, 'b', t_out_orbit, u_prop, 'r--')
ylabel('u [rad]')
xlabel('Orbits')
legend('State Propagation', 'OE Propagation')
grid on

%sgtitle('Short Period State Propagation vs OE Propagation for J2 Effects')

%long period
