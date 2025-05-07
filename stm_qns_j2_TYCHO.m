function [STMs, roe_output] = stm_qns_j2_TYCHO(t, qns_roe_init, SV1_oe_init)
    % Input t is the all the time-steps that we want to propagate the ROE states by
    % qns_roe_init is the initial conditions of the deputy in quasi-nonsingular relative orbital elements
    % SV1_oe_init is the initial conditions of the chief in absolute singular orbital elements
    global J2 Re mu
    a = SV1_oe_init(1);
    ecc = SV1_oe_init(2);
    i_SV1_init = SV1_oe_init(3);
    i = (i_SV1_init);
    w_SV1_init = SV1_oe_init(5);
    w = (w_SV1_init);
    eta = sqrt(1 - ecc^2);
    E = 1 + eta;
    kappa = 3/4*J2*Re^2 * sqrt(mu)/(a^(7/2) * eta^4);
    F = 4 + 3*eta;
    G = 1/eta^2;
    P = 3*(cos(i))^2 - 1;
    Q = 5*(cos(i))^2 - 1;
    R = cos(i);
    S = sin(2*i);
    T = (sin(i))^2;
    U = sin(i);
    V = tan(i/2);
    W = (cos(i/2))^2;
    n = sqrt(mu/a^3);
    omega_dot = kappa * Q;

    exi = ecc*cos(w);
    eyi = ecc*sin(w); 
    
    n_t = length(t);
    roe_output = zeros(n_t, 6);
    STMs = zeros(n_t, 6, 6);
    roe_output(1, :) = qns_roe_init;

    for iter = 1:n_t
        tau = t(iter);
        exf = ecc*cos(omega_dot*tau + w);
        eyf = ecc*sin(omega_dot*tau + w);
        STM = [1,                            0, 0,                                           0,                                            0,                  0;...
            -(3/2*n + 7/2*kappa*E*P)*tau, 1, kappa*exi*F*G*P*tau,                         kappa*eyi*F*G*P*tau,                          -kappa*F*S*tau,     0;...
            7/2*kappa*eyf*Q*tau,          0, cos(omega_dot*tau)-4*kappa*exi*eyf*G*Q*tau, -sin(omega_dot*tau)-4*kappa*eyi*eyf*G*Q*tau, 5*kappa*eyf*S*tau,  0;...
            -7/2*kappa*exf*Q*tau,         0, sin(omega_dot*tau)-4*kappa*exi*exf*G*Q*tau, cos(omega_dot*tau)-4*kappa*eyi*exf*G*Q*tau,  -5*kappa*exf*S*tau, 0;...
            0,                            0, 0,                                           0,                                            1,                  0;...
            7/2*kappa*S*tau,              0, -4*kappa*exi*G*S*tau,                        -4*kappa*eyi*G*S*tau,                         2*kappa*T*tau,      1];
        STMs(iter, :, :) = STM;
        roe_output(iter, :) = STM * qns_roe_init';
    end
end