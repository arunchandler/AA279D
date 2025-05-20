function STM = calc_STM_for_control(T_opt, tf, SV1_oe_init)
    % Use the same globals as your main script
    global J2 Re mu

    % unpack OE (all already in meters/radians)
    a      = SV1_oe_init(1);
    ecc    = SV1_oe_init(2);
    i_rad  = SV1_oe_init(3);

    eta    = sqrt(1 - ecc^2);
    % J2 drift constant
    kappa  = 3/4 * J2 * Re^2 * sqrt(mu) / ( a^(7/2) * eta^4 );

    % inclination‐dependent terms
    P      = 3*cos(i_rad)^2 - 1;
    Q      = 5*cos(i_rad)^2 - 1;
    S      = sin(2*i_rad);
    Tterm  = sin(i_rad)^2;

    n        = sqrt(mu/a^3);
    omega_dot = kappa * Q;
    tau       = tf - T_opt;

    % each row must be exactly 1×6
    STM = [ ...
      1,                                              0,               0,               0,                   0,               0; ...
     -(3/2*n + 7/2 * kappa * (1+eta)*P)*tau,          1,               0,               0,    -kappa * (4+3*eta) * S * tau,     0; ...
      0,                                              0,  cos(omega_dot*tau), -sin(omega_dot*tau),   0,               0; ...
      0,                                              0,  sin(omega_dot*tau),  cos(omega_dot*tau),   0,               0; ...
      0,                                              0,               0,               0,                   1,               0; ...
      7/2 * kappa * S * tau,                          0,               0,               0,    2 * kappa * Tterm * tau,          1 ...
    ];
end
