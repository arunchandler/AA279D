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

    % reduced plant matrix A_c (Eq.6) :contentReference[oaicite:1]{index=1}:contentReference[oaicite:2]{index=2}
    A_c = kappa * [ ...
        0,           0,               0,               0,    0,  1/kappa; ...
        (7/2)*ey*Qp, -(4*ex*ey*G + C)*Qp,  -( (1+4*ey^2*G) - D )*Qp,  5*ey*S, 0,  D - ex; ...
       -(7/2)*ex*Qp,  ( (1+4*ex^2*G) - D )*Qp,  (4*ex*ey*G - C)*Qp, -5*ex*S, 0,  C - ey; ...
        0,           0,               0,               0,    0,  0; ...
        (7/2)*S,    -4*ex*G*S,       -4*ey*G*S,        2*T,   0,  0; ...
        0,           0,               0,               0,    0,  0  ...
    ];

    % reduced input matrix B_c (Eq.7) :contentReference[oaicite:3]{index=3}:contentReference[oaicite:4]{index=4}
    n_c   = sqrt(mu/a^3);
    f_c   = mean2true(oe_c(6), e, tol);         % true anomaly
    denom = a * n_c * (1 + e*cos(f_c));

    B_c = (1/denom) * [ ...
        2/eta,                                 0; ...
        (1+e*cos(f_c)) * (e/eta)*sin(f_c),     0; ...
        (1+e*cos(f_c)) * (e/eta)*cos(f_c),     0; ...
        0,       eta * sin(omega + f_c)/(1+e*cos(f_c)); ...
        0,  eta * ((2+e*cos(f_c))*cos(omega+f_c)+ex)/(1+e*cos(f_c)); ...
        0,  eta * ((2+e*cos(f_c))*sin(omega+f_c)+ey)/(1+e*cos(f_c))  ...
    ];
end