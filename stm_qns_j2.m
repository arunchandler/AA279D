function phi = stm_qns_j2(tau, oe)
%STM_QNS_J2   State transition matrix for quasi-nonsingular ROE with J2
%
%   Φ = stm_qns_j2(τ, oe) returns the 6×6 mapping
%       [δa; δλ; δe_x; δe_y; δi_x; δi_y]₀ → [ ]ₜ
%   under J₂ perturbations, linearized about the chief orbit.
%
%   Inputs:
%     τ   – propagation time [s]
%     oe  – chief classical elements [a; e; i; Ω; ω; M]
%            a = semi-major axis [m]
%            e = eccentricity
%            i = inclination [rad]
%            Ω = RAAN [rad]
%            ω = argument of perigee [rad]
%            M = mean anomaly [rad]
%
%   Output:
%     phi   – 6×6 state transition matrix

global mu J2 Re 

% unpack chief elements
a     = oe(1);
e     = oe(2);
i     = oe(3);
%Ω    = oe(4);  % not needed explicitly here
omega = oe(5);
%M    = oe(6);  % not needed explicitly here

% mean motion
n     = sqrt(mu/a^3);

% eccentricity-dependent parameters
eta     = sqrt(1 - e^2);
k     = 3/4 * J2 * Re^2 * sqrt(mu) / (a^(7/2) * eta^4);
E     = 1 + eta;
F     = 4 + 3*eta;
G     = 1/eta^2;

% inclination-dependent parameters
P     = 3*cos(i)^2 - 1;
Q     = 5*cos(i)^2 - 1;
S     = sin(2*i);
T     = sin(i)^2;

% argument-of-perigee rate
omega_dot = k * Q;
omega_f = omega_dot * tau;
theta         = omega_dot * tau;

% chief e-vector components at t=0
e_xi = e*cos(omega);
e_yi = e*sin(omega);
% rotate them forward by J2-drift
e_xf =  e*cos(omega_f);
e_yf =  e*sin(omega_f);

% build STM
phi = zeros(6,6);

% row 1
phi(1,1) = 1;

% row 2
phi(2,1) = -(1.5*n + 3.5*k*E*P) * tau;
phi(2,2) = 1;
phi(2,3) =  k * e_xi * F * G * P * tau;
phi(2,4) =  k * e_yi * F * G * P * tau;
phi(2,5) = -k * F * S * tau;

% row 3
phi(3,1) =  3.5 * k * e_yf * Q * tau;
phi(3,3) =  cos(theta) - 4*k * e_xi * e_yf * G * Q * tau;
phi(3,4) = -sin(theta) - 4*k * e_yi * e_yf * G * Q * tau;
phi(3,5) =  5 * k * e_yf * S * tau;

% row 4
phi(4,1) = -3.5 * k * e_xf * Q * tau;
phi(4,3) =  sin(theta) + 4*k * e_xi * e_xf * G * Q * tau;
phi(4,4) =  cos(theta) + 4*k * e_yi * e_xf * G * Q * tau;
phi(4,5) = -5 * k * e_xf * S * tau;

% row 5
phi(5,5) = 1;

% row 6
phi(6,1) = 7/2*k *S    *tau;
phi(6,3) = -4*k*e_xi *G* S *tau;
phi(6,4) = -4*k*e_yi *G* S *tau;
phi(6,5) = 2*k     *T    *tau;
phi(6,6) = 1;

end