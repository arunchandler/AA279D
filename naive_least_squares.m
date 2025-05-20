function [delta_v_vals] = naive_least_squares(t_maneuvers, roe_init, roe_final, TSX_oe_init, u_TSX_init, t_0, t_final)
    m = length(t_maneuvers);
    STM_from_init = calc_STM_for_control(t_0, t_final, TSX_oe_init);
    no_control_roe_final = STM_from_init * roe_init';
    total_delta_roe = roe_final' - no_control_roe_final; % Total roe change
    Matrix_block = zeros(6, 3*m);

    for i=1:m
        t_val = t_maneuvers(i);
        STM = calc_STM_for_control(t_val, t_final, TSX_oe_init);
        Gamma = calc_Gamma_for_control(t_val, TSX_oe_init, u_TSX_init);
        Matrix_block(:, 3*(i-1)+1:3*i) = STM * Gamma;
    end
 
    delta_v_vals = reshape(pinv(Matrix_block) * total_delta_roe,m,3)';
    
end