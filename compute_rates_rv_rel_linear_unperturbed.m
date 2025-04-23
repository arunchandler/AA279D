% rates function for numerical integration for propagating rv state
% using the HCW equations

function rv_dot = compute_rates_rv_rel_linear_unperturbed(t, rv)

    global mu;
    global n;
    
    rv_dot = zeros(6,1);
    rv_dot(1:3) = rv(4:6);

    rv_dot(4) = 2*n*rv(5) + 3*n^2*rv(1);
    rv_dot(5) = -2*n*rv(4);
    rv_dot(6) = -n^2*rv(3);

end