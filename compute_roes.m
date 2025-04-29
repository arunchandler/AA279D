function roe = compute_roes(a_ref, i_ref, e_ref, RAAN_ref, u_ref, a, i, e, RAAN, u)
    % ref is chief
    roe = zeros(6,1);
    roe(1) = (a - a_ref) / a_ref;                         % δa
    roe(2) = wrapTo2Pi((u - u_ref) + (RAAN - RAAN_ref) * cos(i_ref)); % δλ
    roe(3) = e * cos(RAAN) - e_ref * cos(RAAN_ref);       % δe_x
    roe(4) = e * sin(RAAN) - e_ref * sin(RAAN_ref);       % δe_y
    roe(5) = i - i_ref;                                   % δi_x
    roe(6) = (RAAN - RAAN_ref) * sin(i_ref);              % δi_y
end