% rates function for numerical integration for relative motion
% unperturbed
% Inputs:
%   state : TDX rv_RTN 6x1, r0 (TSX in ECI) 6x1 
% Outputs: 

function [state_dot] = compute_rates_rv_rel_unperturbed_RTN(t, state)

    global mu;

    r_RTN = state(1:3);
    v_RTN = state(4:6);
    r_TSX = state(7:9); % in ECI
    r_TSX_mag = norm(r_TSX);
    v_TSX = state(10:12);
    v_TSX_mag = norm(v_TSX);
    h = cross(r_TSX, v_TSX);
    h_mag = norm(h);
    
    state_dot = zeros(12,1);
    state_dot(1:3) = v_RTN;
    state_dot(7:9) = v_TSX;
    state_dot(10:12) = -mu .* r_TSX / r_TSX_mag^3; %r0 acceleration in ECI

    theta_dot = h_mag / r_TSX_mag^2;
    r_dot = dot(r_TSX, v_TSX) / r_TSX_mag;
    theta_ddot = -2 * theta_dot * r_dot / r_TSX_mag;
    x = r_RTN(1);
    y = r_RTN(2);
    z = r_RTN(3);
    x_dot = v_RTN(1);
    y_dot = v_RTN(2);

    fact = -mu  / ((r_TSX_mag+x)^2+y^2+z^2)^(3/2);

    R_dd = 2*theta_dot*y_dot + theta_ddot*y + theta_dot^2*x + (r_TSX_mag+x)*fact + mu/r_TSX_mag^2; 
    T_dd = -2*theta_dot*x_dot - theta_ddot*x + theta_dot^2*y + y*fact;
    N_dd = z*fact;

    state_dot(4:6) = [R_dd;T_dd;N_dd];

end