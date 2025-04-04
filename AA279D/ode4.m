function [t_out, y_out] = ode4(odefun, tspan, y0, h)
    % Fixed-step 4th-order Runge-Kutta solver
    t_out = tspan(1):h:tspan(2);
    n_steps = length(t_out);
    y_out = zeros(n_steps, length(y0));
    y_out(1, :) = y0;
    y = y0;
    
    for i = 1:n_steps-1
        t = t_out(i);
        k1 = h * odefun(t, y);
        k2 = h * odefun(t + 0.5*h, y + 0.5*k1);
        k3 = h * odefun(t + 0.5*h, y + 0.5*k2);
        k4 = h * odefun(t + h, y + k3);
        
        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        y_out(i+1, :) = y;
    end
end