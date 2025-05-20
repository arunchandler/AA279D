function [A_c, B_c] = plant_qns(oe_c)
%PLANT_QNS QNS plant & input matrices (no drag)
%   [A_c,B_c] = plant_reduced_qns(oe_c)
%   oe_c = [a; e; i; RAAN; omega; M]  (chief mean elements)
%
%   A_c (7×7) and B_c (7×3) implement
%     [δα̇; δä] = A_c [δα; δȧ] + [B_c; 0(1×2)]·u
%   for δα = [δa; δlambda; δex; δey; δix; δiy] per Steindorf Eq.(A.1 & A.7).

    global mu J2 Re tol

    % unpack
    a      = oe_c(1);
    e      = oe_c(2);
    i      = oe_c(3);
    omega  = oe_c(5);

    % common substitutions
    eta    = sqrt(1 - e^2);
    gamma  = 3/4 * J2 * Re^2 * sqrt(mu);
    kappa  = gamma / ( a^(7/2) * eta^4 );

    ex     = e * cos(omega);
    ey     = e * sin(omega);
    C      = sin(omega);
    D      = cos(omega);
    E      = 1 + eta;
    F      = 4 + 3*eta;
    G      = 1/eta^2;
    Pp     = 3*cos(i)^2 - 1;
    Qp     = 5*cos(i)^2 - 1;
    S      = sin(2*i);
    T      = sin(i)^2;

    n_c   = sqrt(mu/a^3);
    f_c   = mean2true(oe_c(6), e, tol);         % true anomaly
    th = f_c + omega;
    term1 = (1 + e*cos(f_c));
    term2 = 2+e*cos(f_c);

    % --- build Ā (7×7) ---
    Abar = zeros(7,7);

    % row 2
    Abar(2,:) = [
        -7/2 * E * Pp, ...
         0, ...
         ex * G * F * Pp, ...
         ey * G * F * Pp, ...
        -F * S, ...
         0, ...
         0 ];

    % row 3
    Abar(3,:) = [
         7/2 * ey * Qp, ...
         0, ...
        -4 * ex * ey * G * Qp, ...
        -(1 + 4*ex^2 * G)*Qp, ...
         5 * ey * S, ...
         0, ...
         0 ];

    % row 4
    Abar(4,:) = [
        -7/2 * ex * Qp, ...
         0, ...
         (1 + 4*ex^2 * G)*Qp, ...
         4 * ex * ey * G * Qp, ...
        -5 * ex * S, ...
         0, ...
         0 ];

    % row 6
    Abar(6,:) = [
         7/2 * S, ...
         0, ...
        -4 * ex * G * S, ...
        -4 * ey * G * S, ...
         2 * T, ...
         0, ...
         0 ];

    % the other rows (1,5,7) remain zero
    A_c = kappa * Abar;

    % build 6×3 B̄
    Bbar = zeros(6,3);
    Bbar(1,1) = 2*e*sin(f_c)      / eta;      % row 1
    Bbar(1,2) = 2*term1           / eta;

    Bbar(2,1) = -2*eta^2           / term1;    % row 2

    Bbar(3,1) =  eta * sin(th);                 % row 3
    Bbar(3,2) =  eta * ( term2*cos(th) + ex ) / term1;
    Bbar(3,3) =  eta * ey * sin(th)          / ( tan(i)*term1 );

    Bbar(4,1) = -eta * cos(th);                 % row 4
    Bbar(4,2) =  eta * ( term2*sin(th) + ey ) / term1;
    Bbar(4,3) = -eta * ex * sin(th)          / ( tan(i)*term1 );

    Bbar(5,3) =  eta * cos(th) / term1;         % row 5
    Bbar(6,3) =  eta * sin(th) / term1;         % row 6

    % append zero‐row and scale by 1/(a n_c)
    B_c = (1/(a*n_c)) * [ Bbar; zeros(1,3) ];

end