function Gamma = calc_Gamma_for_control(t_man, SV1_oe_init, u_SV1_init)
    global J2 Re mu

    a     = SV1_oe_init(1);
    ecc   = SV1_oe_init(2);
    i_rad = SV1_oe_init(3);

    eta   = sqrt(1 - ecc^2);
    kappa = 3/4 * J2 * Re^2 * sqrt(mu) / ( a^(7/2) * eta^4 );

    % mean motion
    n     = sqrt(mu/a^3);

    % argument of latitude in radians (u_SV1_init must already be rad)
    u = t_man*(n + kappa*(eta*(3*cos(i_rad)^2-1) + (5*cos(i_rad)^2-1))) ...
        + u_SV1_init;

    % each row = 1×3, so Gamma is 6×3
    Gamma = 1/(n*a) * [ ...
       0,       2,      0; ...
      -2,       0,      0; ...
       sin(u),  2*cos(u), 0; ...
      -cos(u),  2*sin(u), 0; ...
       0,       0,      cos(u); ...
       0,       0,      sin(u) ...
    ];
end
