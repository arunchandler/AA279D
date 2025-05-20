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
tend = s_d*0.5;
t_vec = [tstart:tint:tend];
num_points = floor((tend-tstart)/tint) + 1;


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

n = sqrt(mu/a_TSX_2^3);
T = 2*pi/n;

t_out_orbit = t_vec / T;

%% Part b) numerical integration of nonlinear equations of relative motion

[init_TDX_rv_RTN,~] = eci2rtn(init_TSX_rv_ECI, init_TDX_rv_ECI);
r0_init = init_TSX_rv_ECI;

init_state = [init_TDX_rv_RTN; r0_init];

[t_out, TDX_rv_out_unperturbed_nonlin] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [tstart, tend]', init_state, tint);

figure;
%sgtitle('Relative Position & Velocity');
subplot(2,3,1)
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,1));
xlabel('Orbits');
ylabel('R Position [m]');
grid on;

subplot(2,3,2)
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,2));
xlabel('Orbits');
ylabel('T Position [m]');
grid on;

subplot(2,3,3)
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,3));
xlabel('Orbits');
ylabel('N Position [m]');
grid on;

subplot(2, 3, 4)
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:, 4));
xlabel('Orbits');
ylabel('R Velocity [m/s]');
grid on;

subplot(2, 3, 5)
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:, 5));
xlabel('Orbits');
ylabel('T Velocity [m/s]');
grid on;

subplot(2, 3, 6)
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:, 6));
xlabel('Orbits');
ylabel('N Velocity [m/s]');
grid on;

figure;
subplot(2,2,1);
plot(TDX_rv_out_unperturbed_nonlin(:,2), TDX_rv_out_unperturbed_nonlin(:,1));
grid on;
axis equal;
xlabel('T Position [m]');
ylabel('R Position [m]');

subplot(2,2,2);
plot(TDX_rv_out_unperturbed_nonlin(:,3), TDX_rv_out_unperturbed_nonlin(:,1));
grid on;
axis equal;
xlabel('N Position [m]');
ylabel('R Position [m]');

subplot(2,2,3);
plot(TDX_rv_out_unperturbed_nonlin(:,3), TDX_rv_out_unperturbed_nonlin(:,2));
grid on;
axis equal;
xlabel('N Position [m]');
ylabel('T Position [m]');

subplot(2,2,4)
plot3(TDX_rv_out_unperturbed_nonlin(:,1), TDX_rv_out_unperturbed_nonlin(:,2), TDX_rv_out_unperturbed_nonlin(:,3));
hold on;
scatter3(0, 0, 0, 100, 'Marker', '.');
grid on;
axis equal;
xlabel('R Position [m]');
ylabel('T Position [m]');
zlabel('N Position [m]');
%legend({'One','Two','Three','Four'},'FontSize',14)
legend({'TDX Position', 'TSX Position'},'FontSize',8)
%title('3D Relative Orbit of TDX in TSX''s RTN Frame');

%% Part c) compute relative orbit with fundamental diff eq

[t_out, TSX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TSX_rv_ECI, tint);
[t_out, TDX_rv_out_unperturbed] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TDX_rv_ECI, tint);

TDX_rv_RTN = zeros(num_points, 6);

for idx = 1:num_points
    TSX_rv = TSX_rv_out_unperturbed(idx, :)';
    TDX_rv = TDX_rv_out_unperturbed(idx, :)';
    [TDX_rv_RTN(idx, :), ~] = eci2rtn(TSX_rv, TDX_rv);
end

figure;
%sgtitle('Relative Position & Velocity Comparison');
subplot(2,3,1)
plot(t_out_orbit, TDX_rv_RTN(:,1));
hold on;
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,1));
xlabel('Orbits');
ylabel('R Position [m]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,2)
plot(t_out_orbit, TDX_rv_RTN(:,2));
hold on;
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,2));
xlabel('Orbits');
ylabel('T Position [m]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,3)
plot(t_out_orbit, TDX_rv_RTN(:,3));
hold on;
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,3));
xlabel('Orbits');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
ylabel('N Position [m]');
grid on;

subplot(2,3,4)
plot(t_out_orbit, TDX_rv_RTN(:,4));
hold on;
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,4));
xlabel('Orbits');
ylabel('R Velocity [m/s]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,5)
plot(t_out_orbit, TDX_rv_RTN(:,5));
hold on;
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,5));
xlabel('Orbits');
ylabel('T Velocity [m/s]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

subplot(2,3,6)
plot(t_out_orbit, TDX_rv_RTN(:,6));
hold on;
plot(t_out_orbit, TDX_rv_out_unperturbed_nonlin(:,6));
xlabel('Orbits');
ylabel('N Velocity [m/s]');
legend('Fundamental Differential Eq.', 'Relative Nonlinear Eq.');
grid on;

figure;

subplot(2,2,1);
plot(TDX_rv_RTN(:,2), TDX_rv_RTN(:,1));
hold on;
plot(TDX_rv_out_unperturbed_nonlin(:,2), TDX_rv_out_unperturbed_nonlin(:,1));
grid on;
axis equal;
xlabel('T Position [m]');
ylabel('R Position [m]');

subplot(2,2,2);
plot(TDX_rv_RTN(:,3), TDX_rv_RTN(:,1));
hold on;
plot(TDX_rv_out_unperturbed_nonlin(:,3), TDX_rv_out_unperturbed_nonlin(:,1));
grid on;
axis equal;
xlabel('N Position [m]');
ylabel('R Position [m]');

subplot(2,2,3);
plot(TDX_rv_RTN(:,3), TDX_rv_RTN(:,2));
hold on;
plot(TDX_rv_out_unperturbed_nonlin(:,3), TDX_rv_out_unperturbed_nonlin(:,2));
grid on;
axis equal;
xlabel('N Position [m]');
ylabel('T Position [m]');

subplot(2,2,4);
plot3(TDX_rv_RTN(:,1), TDX_rv_RTN(:,2), TDX_rv_RTN(:,3));
hold on;
plot3(TDX_rv_out_unperturbed_nonlin(:,1), TDX_rv_out_unperturbed_nonlin(:,2), TDX_rv_out_unperturbed_nonlin(:,3));
scatter3(0, 0, 0, 100, 'Marker', '.');
grid on;
xlabel('R Position [m]');
ylabel('T Position [m]');
zlabel('N Position [m]');
axis equal;
legend({'TDX Position (fundamental diff eq.)', 'TDX Position (relative nonlin. eq.)', 'TSX Position'}, 'FontSize',8)
%sgtitle('3D Relative Orbit of TDX in TSX''s RTN Frame');

%% Part d) compare b and c

err_RTN = TDX_rv_out_unperturbed_nonlin(:, 1:6) - TDX_rv_RTN;

figure;
subplot(2,3,1)
plot(t_out_orbit, err_RTN(:,1));
xlabel('Orbits');
ylabel('Error in R Position [m]');
grid on;

subplot(2,3,2)
plot(t_out_orbit, err_RTN(:,2));
xlabel('Orbits');
ylabel('Error in T Position [m]');
grid on;

subplot(2,3,3)
plot(t_out_orbit, err_RTN(:,3));
xlabel('Orbits');
ylabel('Error in N Position [m]');
grid on;

subplot(2,3,4)
plot(t_out_orbit, err_RTN(:,4));
xlabel('Orbits');
ylabel('Error in R Velocity [m/s]');
grid on;

subplot(2,3,5)
plot(t_out_orbit, err_RTN(:,5));
xlabel('Orbits');
ylabel('Error in T Velocity [m/s]');
grid on;

subplot(2,3,6)
plot(t_out_orbit, err_RTN(:,6));
xlabel('Orbits');
ylabel('Error in N Velocity [m/s]');
grid on;

%sgtitle('Error Comparison: Unperturbed Nonlinear Eq. vs. Differential ECI');

% part d) second piece with a non-zero difference in semi-major axis 
a_TSX_2_d = 6886536.686;  % m
i_TSX_2_d = deg2rad(97.4453);
e_TSX_2_d = 0.0001264;
RAAN_TSX_2_d = deg2rad(351.0108);
omega_TSX_2_d = deg2rad(101.2452);
M_TSX_2_d = deg2rad(11.6520);
nu_TSX_2_d = mean2true(M_TSX_2, e_TSX_2,tol);
u_TSX_2_d = nu_TSX_2 + omega_TSX_2;

a_TDX_2_d =  6881866.120352;  % changed back to original semi-major axis so there is an offset
i_TDX_2_d = deg2rad(97.4454);
e_TDX_2_d = 0.0001269;
RAAN_TDX_2_d = deg2rad(351.0106);
omega_TDX_2_d = deg2rad(100.5043);
M_TDX_2_d = deg2rad(12.3926);

init_TSX_koe_d = [a_TSX_2_d, e_TSX_2_d, i_TSX_2_d, RAAN_TSX_2_d, omega_TSX_2_d, M_TSX_2_d];
init_TDX_koe_d = [a_TDX_2_d, e_TDX_2_d, i_TDX_2_d, RAAN_TDX_2_d, omega_TDX_2_d, M_TDX_2_d];

init_TSX_rv_ECI_d = oe2rv(init_TSX_koe_d, mu);
init_TDX_rv_ECI_d = oe2rv(init_TDX_koe_d, mu);

[init_TDX_rv_RTN_d,~] = eci2rtn(init_TSX_rv_ECI_d, init_TDX_rv_ECI_d);
r0_init_d = init_TSX_rv_ECI_d;

init_state_d = [init_TDX_rv_RTN_d; r0_init_d];
% relative propagation
[t_out, TDX_rv_out_unperturbed_nonlin_d] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [tstart, tend]', init_state_d, tint);

%fundamental eq propagation

[t_out, TSX_rv_out_unperturbed_d] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TSX_rv_ECI_d, tint);
[t_out, TDX_rv_out_unperturbed_d] = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', init_TDX_rv_ECI_d, tint);

TDX_rv_RTN_d = zeros(num_points, 6);

for idx = 1:num_points
    TSX_rv_d = TSX_rv_out_unperturbed_d(idx, :)';
    TDX_rv_d = TDX_rv_out_unperturbed_d(idx, :)';
    [TDX_rv_RTN_d(idx, :), ~] = eci2rtn(TSX_rv_d, TDX_rv_d);
end

err_RTN_d = TDX_rv_out_unperturbed_nonlin_d(:, 1:6) - TDX_rv_RTN_d;

figure;
subplot(2,3,1)
plot(t_out_orbit, err_RTN_d(:,1));
xlabel('Orbits');
ylabel('Error in R Position [m]');
grid on;

subplot(2,3,2)
plot(t_out_orbit, err_RTN_d(:,2));
xlabel('Orbits');
ylabel('Error in T Position [m]');
grid on;

subplot(2,3,3)
plot(t_out_orbit, err_RTN_d(:,3));
xlabel('Orbits');
ylabel('Error in N Position [m]');
grid on;

subplot(2,3,4)
plot(t_out_orbit, err_RTN_d(:,4));
xlabel('Orbits');
ylabel('Error in R Velocity [m/s]');
grid on;

subplot(2,3,5)
plot(t_out_orbit, err_RTN_d(:,5));
xlabel('Orbits');
ylabel('Error in T Velocity [m/s]');
grid on;

subplot(2,3,6)
plot(t_out_orbit, err_RTN_d(:,6));
xlabel('Orbits');
ylabel('Error in N Velocity [m/s]');
grid on;

%sgtitle('Error Comparison: Unperturbed Nonlinear Eq. vs. Differential ECI with Semi-major Axis Offset');


%% Part e) most fuel-efficient impulse maneuver by deputy - try integrating analytically the GVEs over a dv
delta_a = a_TSX_2_d - a_TDX_2_d;

d_v = (sqrt(mu)/(2*a_TSX_2_d^(3/2)))*delta_a;

%% Part f) Apply maneuver by repeating b and adding discontinuity in the inertial velocity of deputy

M0 = M_TDX_2_d;
n_TDX = sqrt(mu/(a_TDX_2_d^3));
M_arr = wrapTo2Pi(M0 + n_TDX*t_out);

%loop through M_arr to find periapsis points
t_out = t_out(:);

t_m_candidate = 1e4;

indices_before = find(t_out <= t_m_candidate);

[~, local_idx] = min(M_arr(indices_before));  
idx_perigee = indices_before(local_idx);
t_perigee = t_out(idx_perigee);

%T_orbit = 2*pi*sqrt(a_TSX_2_d^3/mu);

n_orbits = 3;  
t_m_new = t_perigee + n_orbits * T;

fprintf('Chosen perigee time: %g seconds\n', t_perigee);
fprintf('New maneuver time (after %d orbits): %g seconds\n', n_orbits, t_m_new);

[t_out_pre, state_out_pre] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [tstart, t_m_new]', init_state_d, tint);
t_out_pre = t_out_pre(:);
state_at_tm = state_out_pre(end, :)';


impulse_vector = [0; d_v; 0];
state_at_tm(4:6) = state_at_tm(4:6) + impulse_vector;


[t_out_post, state_out_post] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [t_m_new, tend]', state_at_tm, tint);
t_out_post = t_out_post(:);


t_out_d_corrected = [t_out_pre; t_out_post(2:end)];
TDX_rv_out_unperturbed_nonlin_d_corrected = [state_out_pre; state_out_post(2:end,:)];

t_out_orbit_no_maneuver = t_out / T;
t_out_orbit_corrected = t_out_d_corrected / T;

figure;
%sgtitle('Relative Position & Velocity Comparison (_d Case)');

% R Position
subplot(2,3,1)
plot(t_out_orbit_no_maneuver, TDX_rv_out_unperturbed_nonlin_d(:,1), 'b-', 'LineWidth',2);
hold on;
plot(t_out_orbit_corrected, TDX_rv_out_unperturbed_nonlin_d_corrected(:,1), 'r-', 'LineWidth',2);
xlabel('Orbits');
ylabel('R Position [m]');
legend('No Maneuver', 'Corrected','Location','Best');
grid on;

% T Position
subplot(2,3,2)
plot(t_out_orbit_no_maneuver, TDX_rv_out_unperturbed_nonlin_d(:,2), 'b-', 'LineWidth',2);
hold on;
plot(t_out_orbit_corrected, TDX_rv_out_unperturbed_nonlin_d_corrected(:,2), 'r-', 'LineWidth',2);
xlabel('Orbits');
ylabel('T Position [m]');
legend('No Maneuver', 'Corrected','Location','Best');
grid on;

% N Position
subplot(2,3,3)
plot(t_out_orbit_no_maneuver, TDX_rv_out_unperturbed_nonlin_d(:,3), 'b-', 'LineWidth',2);
hold on;
plot(t_out_orbit_corrected, TDX_rv_out_unperturbed_nonlin_d_corrected(:,3), 'r-', 'LineWidth',2);
xlabel('Orbits');
ylabel('N Position [m]');
legend('No Maneuver', 'Corrected','Location','Best');
grid on;

% R Velocity
subplot(2,3,4)
plot(t_out_orbit_no_maneuver, TDX_rv_out_unperturbed_nonlin_d(:,4), 'b-', 'LineWidth',2);
hold on;
plot(t_out_orbit_corrected, TDX_rv_out_unperturbed_nonlin_d_corrected(:,4), 'r-', 'LineWidth',2);
xlabel('Orbits');
ylabel('R Velocity [m/s]');
legend('No Maneuver', 'Corrected','Location','Best');
grid on;

% T Velocity
subplot(2,3,5)
plot(t_out_orbit_no_maneuver, TDX_rv_out_unperturbed_nonlin_d(:,5), 'b-', 'LineWidth',2);
hold on;
plot(t_out_orbit_corrected, TDX_rv_out_unperturbed_nonlin_d_corrected(:,5), 'r-', 'LineWidth',2);
xlabel('Orbits');
ylabel('T Velocity [m/s]');
legend('No Maneuver', 'Corrected','Location','Best');
grid on;

% N Velocity
subplot(2,3,6)
plot(t_out_orbit_no_maneuver, TDX_rv_out_unperturbed_nonlin_d(:,6), 'b-', 'LineWidth',2);
hold on;
plot(t_out_orbit_corrected, TDX_rv_out_unperturbed_nonlin_d_corrected(:,6), 'r-', 'LineWidth',2);
xlabel('Orbits');
ylabel('N Velocity [m/s]');
legend('No Maneuver', 'Corrected','Location','Best');
grid on;


burn_orbit = t_m_new / T;

window = 0.2;  
xlim_min = burn_orbit - window;
xlim_max = burn_orbit + window;


figure;

plot3(TDX_rv_out_unperturbed_nonlin_d(:,1), TDX_rv_out_unperturbed_nonlin_d(:,2), TDX_rv_out_unperturbed_nonlin_d(:,3), 'b-', 'LineWidth', 2);
hold on;
plot3(TDX_rv_out_unperturbed_nonlin_d_corrected(:,1), TDX_rv_out_unperturbed_nonlin_d_corrected(:,2), TDX_rv_out_unperturbed_nonlin_d_corrected(:,3), 'r-', 'LineWidth', 2);
scatter3(0, 0, 0, 100, 'k.', 'MarkerEdgeColor','k');  % Chief's position at origin
grid on;
%axis equal;
xlabel('R Position [m]');
ylabel('T Position [m]');
zlabel('N Position [m]');
%title('3D Relative Orbit');
legend({'Fundamental diff eq.', 'Corrected Simulation', 'Chief (TSX)'}, 'FontSize', 8);