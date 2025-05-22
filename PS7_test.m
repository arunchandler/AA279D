clc; clear; close all;
addpath('mean_osc');
format long g;

%-----------------------
% Propagate non‑linear + STM in one go,
% then reconstruct linearized perturbed traj.
%-----------------------

global tol Re J2 mu s_d sigma
tol = 1e-10;
Re  = 6378137;          % m
J2  = 1.082626e-3;
mu  = 3.986004418e14;   % m^3/s^2
s_d = 86400;            % s/day
sigma = 1e-9;

%--- initial chief orbit in rv form ---
a_init    = 6886536.686;  % m
e_init    = 0.0001264;
i_init    = deg2rad(97.4453);
RAAN_init = deg2rad(351.0108);
omega_init= deg2rad(101.2452);
M_init    = deg2rad(11.6520);
oe_init   = [a_init,e_init,i_init,RAAN_init,omega_init,M_init];
rv_init   = oe2rv(oe_init,mu);

% timing parameters
tstart     = 0.0;
n          = sqrt(mu/a_init^3);
T          = 2*pi/n;
n_orbit    = 15;
tend       = n_orbit*T;
num_points = 10000;
dt         = (tend - tstart)/(num_points-1);
t_grid     = linspace(tstart, tend, num_points).';
t_orbit    = t_grid / T;
orbit_num = floor(t_orbit) + 1;

delta0     = [0; 0; 0; 0; 0; 0];  % 10 m radial error, zero else

%%── Stack initial state + STM into one vector X0 ────────────────────
Phi0 = eye(6);
X0   = [rv_init; Phi0(:)];        % 6 + 36 = 42 elements

%%── Integrate the augmented system once ─────────────────────────────
%   - first 6 components           = true nonlinear [r;v]
%   - next 36 components (Φ as vec) = state‑transition matrix
[t_out, rv_nonlin] = ode4(@compute_rates_rv_perturbed, [tstart,tend]', rv_init, dt);
[t_out, X] = ode4(@compute_rates_aug, [tstart,tend]', X0, dt);

%%── Unpack results ─────────────────────────────────────────────────
rv_nom      = X(:,1:6);           % nominal nonlinear trajectory
Phi_flat    = X(:,7:end);         % each row is 1×36 = vec(Φ)

% reconstruct linearized trajectory: x_lin(t) = x_nom(t) + Φ(t)*δ₀
rv_lin = zeros(size(rv_nom));
for k = 1:size(X,1)
    Phi_k       = reshape(Phi_flat(k,:),6,6);
    rv_lin(k,:) = (rv_nom(k,:).' + Phi_k*delta0).';
end


%Plot in 3D
figure;
plot3(rv_nonlin(:,1), rv_nonlin(:,2), rv_nonlin(:,3), 'b', 'DisplayName','Nonlinear');
hold on;
plot3(rv_lin(:,1)   , rv_lin(:,2)   , rv_lin(:,3)   , 'r--','DisplayName','Linearized w Noise');
axis equal;
grid on;
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('3D Trajectory: Nonlinear vs Linearized w Noise');
legend('Location','best');

% compute error in the three position components
pos_err = rv_lin(:,1:3) - rv_nonlin(:,1:3);

% plot errors in one figure with 3 subplots
figure;
% X‐error
subplot(3,1,1);
plot(t_orbit, pos_err(:,1));
xlabel('Orbit Number');
ylabel('x_{lin} - x_{nonlin} (m)');
title('Position Error in X');
grid on;

% Y‐error
subplot(3,1,2);
plot(t_orbit, pos_err(:,2));
xlabel('Orbit Number');
ylabel('y_{lin} - y_{nonlin} (m)');
title('Position Error in Y');
grid on;

% Z‐error
subplot(3,1,3);
plot(t_orbit, pos_err(:,3));
xlabel('Orbit Number');
ylabel('z_{lin} - z_{nonlin} (m)');
title('Position Error in Z');
grid on;