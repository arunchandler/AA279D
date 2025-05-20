function roe = compute_roes(oe_ref, oe)
%compute_roes  Quasi‑nonsingular relative orbital elements (ROEs)
%
%   roe = compute_roes(oe_ref, oe) returns the 6×1 vector
%     [δa; δλ; δeₓ; δeᵧ; δiₓ; δiᵧ] 
%   between a deputy (oe) and a chief (oe_ref), where each is
%     [a, e, i, RAAN, ω, M].
%
%   Requires mean2true(M, e, tol) and global tol defined elsewhere.

    % bring in your mean2true tolerance
    global tol

    % unpack chief
    a_ref     = oe_ref(1);
    e_ref     = oe_ref(2);
    i_ref     = wrapTo2Pi(oe_ref(3));
    RAAN_ref  = wrapTo2Pi(oe_ref(4));
    omega_ref = wrapTo2Pi(oe_ref(5));
    M_ref     = wrapTo2Pi(oe_ref(6));
    u_ref     = wrapTo2Pi(omega_ref + M_ref);

    % unpack deputy
    a     = oe(1);
    e     = oe(2);
    i     = wrapTo2Pi(oe(3));
    RAAN  = wrapTo2Pi(oe(4));
    omega = wrapTo2Pi(oe(5));
    M     = wrapTo2Pi(oe(6));
    u     = wrapTo2Pi(omega + M);

    % relative ROEs
    delta_a      = (a - a_ref) / a_ref;
    delta_lambda = (u - u_ref) + (RAAN - RAAN_ref)*cos(i_ref);
    delta_ex     = e*cos(omega)     - e_ref*cos(omega_ref);
    delta_ey     = e*sin(omega)     - e_ref*sin(omega_ref);
    delta_ix     = i - i_ref;
    delta_iy     = (RAAN - RAAN_ref)*sin(i_ref);

    roe = [ delta_a;
            delta_lambda;
            delta_ex;
            delta_ey;
            delta_ix;
            delta_iy ];
end
