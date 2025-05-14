function x_dot = augmented_dynamics(t, x, delta_nom, k, N_ip, N_oop, u_max)
    global tol;
    % x = [delta(6×1)]
    delta = x(1:6);

    % 2) build plant
    [A_c, B_c] = plant_reduced_qns(oe);

    % 3) Lyapunov shaping
    delta5   = delta(1:5);
    Delta5   = delta5 - delta_nom(1:5);
    phi_ip   = atan2(delta5(3), delta5(2));
    phi_oop  = atan2(delta5(5), delta5(4));
    phi_full = oe(5) + mean2true(oe(6), oe(2), tol);

    Jp = phi_full - phi_ip;
    Hp = phi_full - phi_oop;
    P5 = (1/k) * diag([cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Jp)^N_ip; cos(Hp)^N_oop; cos(Hp)^N_oop]);

    % 4) compute control
    A5 = A_c(1:5,1:5);
    B5 = B_c(1:5,:);
    u  = - pinv(B5) * (A5*delta5 + P5*Delta5);
    % optional: saturate u here

    % 5) propagate delta under thrust
    delta_dot = A_c*delta + B_c*u;

    % pack
    x_dot = delta_dot;
end