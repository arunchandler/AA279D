function delta_v_vals = naive_least_squares( ...
        t_maneuvers, roe_init, roe_final, ...
        SV1_oe_init, u_SV1_init, t_0, t_final)

    m = numel(t_maneuvers);

    % STM from t0→tf, no control
    STM_0f = calc_STM_for_control(t_0, t_final, SV1_oe_init);
    roe_no_ctrl = STM_0f * roe_init;          % unitless

    % total required ΔROE
    total_delta_roe = roe_final - roe_no_ctrl;

    % build the big 6×(3m) matrix
    Matrix_block = zeros(6, 3*m);
    for i = 1:m
        t_abs   = t_maneuvers(i);
        dt_i    = t_abs - t_0; 
        STM_if = calc_STM_for_control(t_abs, t_final, SV1_oe_init);
        Gamma_i = calc_Gamma_for_control(t_abs, SV1_oe_init, u_SV1_init);
        Matrix_block(:, 3*(i-1)+1 : 3*i) = STM_if * Gamma_i;
    end

    % solve for the 3m×1 vector of Δv's
    dv_vec = pinv(Matrix_block) * total_delta_roe;  % [Δv_R1; Δv_T1; Δv_N1; …]

    % reshape into m×3: each row = [Δv_R, Δv_T, Δv_N]
    delta_v_vals = reshape(dv_vec, 3, m)';
end
