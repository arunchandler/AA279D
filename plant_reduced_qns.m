function [A_c, B_c] = plant_reduced_qns(oe_c)
%PLANT_REDUCED_QNS  Reduced QNS plant & input matrices (no drag)
%   [A_c,B_c] = plant_reduced_qns(oe_c)
%   oe_c = [a; e; i; RAAN; omega; M]  (chief mean elements)
%
%   A_c (6×6) and B_c (6×2) implement
%     [δα̇; δä] = A_c [δα; δȧ] + [B_c; 0(1×2)]·u
%   for δα = [δa; δex; δey; δix; δiy] per Steindorf Eq.(5–7).

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

    A_c = kappa * [ ...
        0,           0,               0,               0,    0,  0; ...
        (7/2)*ey*Qp, -(4*ex*ey*G)*Qp,  -((1+4*ey^2*G))*Qp,  5*ey*S, 0,  0; ...
       -(7/2)*ex*Qp,  ((1+4*ex^2*G))*Qp,  (4*ex*ey*G)*Qp, -5*ex*S, 0,  0; ...
        0,           0,               0,               0,    0,  0; ...
        (7/2)*S,    -4*ex*G*S,       -4*ey*G*S,        2*T,   0,  0; ...
        0,           0,               0,               0,    0,  0  ...
    ];

    n_c   = sqrt(mu/a^3);
    f_c   = mean2true(oe_c(6), e, tol);         % true anomaly
    th = f_c + omega;
    term1 = (1 + e*cos(f_c));
    term2 = 2+e*cos(f_c);

    B_c = (1/(a*n_c)) * [ ...
        (2/eta)*term1,                                 0; ...
        eta*(term2*cos(th)+ex)/term1,     eta*ey*sin(th)/(tan(i)*term1); ...
        eta*(term2*sin(th)+ey)/term1,     -eta*ex*sin(th)/(tan(i)*term1); ...
        0,       eta * cos(th)/term1; ...
        0,  eta * sin(th)/term1; ...
        0,  0  ...
    ];

end