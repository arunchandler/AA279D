function [A_kep, A_J2] = plant_matrices_qns(oe_c)
%PLANT_MATRICES_QNS Compute Keplerian and J2 plant matrices for QNS ROEs
%   [A_kep, A_J2] = plant_matrices_qns(oe_c)
%   oe_c = [a, e, i, RAAN, omega, M]

    global mu J2 Re

    % Compute Keplerian and first-order J2 plant matrices
    A_kep = A_kep_qns(oe_c);
    A_J2  = A_J2_qns0(oe_c);
end

function A = A_kep_qns(oe)
% Keplerian plant matrix for QNS ROEs
a    = oe(1);
n    = sqrt(mu / a^3);
A    = zeros(6);
A(2,1) = -1.5 * n;
end

function A = A_J2_qns0(oe)
% First-order J2 perturbation plant matrix in QNS basis
a     = oe(1);
e     = oe(2);
i     = oe(3);
eta   = sqrt(1 - e^2);
kappa = (3 * J2 * Re^2 * sqrt(mu)) / (4 * a^(7/2) * eta^4);

% Auxiliary substitutions (E, F, G, P, Q, S, T)
E = 1 + eta;
G = 1/eta^2;
F = 4 + 3*eta;
P = 3*cos(i)^2 - 1;
Q = 5*cos(i)^2 - 1;
S = sin(2*i);
T = sin(i)^2;

A = zeros(6);
% In-plane rows
A(2,1) = -7/2 * kappa * E * P;
A(2,3) =     kappa * e * F * G * P;
A(2,5) =    -kappa * F * S;
% In-track rows
A(4,1) = -7/2 * kappa * e * Q;
A(4,3) =     kappa * 4 * e^2 * G * Q;
A(4,5) =    -kappa * 5 * e * S;
% Cross-track rows
A(6,1) =  7/2 * kappa * S;
A(6,3) =    -4    * kappa * e * G * S;
A(6,5) =     2    * kappa * T;
end
