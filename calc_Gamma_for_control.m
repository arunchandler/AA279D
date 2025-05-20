function Gamma = calc_Gamma_for_control(t_man, SV1_oe_init, u_SV1_init)
    % Some relevant parameters that are useful
    global J2 Re mu
    a = SV1_oe_init(1);
    ecc = SV1_oe_init(2);
    i_SV1_init = SV1_oe_init(3);
    i = i_SV1_init;
    eta = sqrt(1 - ecc^2);
    kappa = 3/4*J2*Re^2 * sqrt(mu)/(a^(7/2) * eta^4);
    P = 3*(cos(i))^2 - 1;
    Q = 5*(cos(i))^2 - 1;
    n = sqrt(mu/a^3);

    u = t_man*(n + kappa*(eta*P + Q)) + u_SV1_init; % TODO: Need to check this conversion but I think it's ok
    Gamma = 1/(n*a) * [0, 2, 0;...
                       -2, 0, 0;...
                       sin(u), 2*cos(u), 0;...
                       -cos(u), 2*sin(u), 0;...
                       0, 0, cos(u);...
                       0, 0, sin(u);];
    
end