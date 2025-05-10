function [oe_dot] = compute_rates_GVE_J2(~, oe)
    % oe = [a; e; i; Omega; omega; M]
    % returns oe_dot = [da; de; di; dOmega; domega; dM]
    
    global mu J2 Re

    % unpack
    a     = oe(1);
    e     = oe(2);
    i     = oe(3);
    % Omega = oe(4);  % not used in rates
    % omega = oe(5);  % not used in rates
    % M     = oe(6);  % not used in rates

    % mean motion
    n = sqrt(mu / a^3);

    % common factor (Re/(a(1-e^2)))^2
    fac = (Re / ( a*(1 - e^2) ))^2;
    cosi = cos(i);

    % secular J2 rates
    dOmega = -3/2 * n * J2 * fac * cosi;
    domega =  3/4 * n * J2 * fac * (5*cosi^2 - 1);
    dM     =  n + 3/4 * n * J2 * fac * sqrt(1 - e^2) * (3*cosi^2 - 1);

    % no secular change in a, e, i
    da = 0;
    de = 0;
    di = 0;

    % pack
    oe_dot = [ da;
               de;
               di;
               dOmega;
               domega;
               dM ];
end
