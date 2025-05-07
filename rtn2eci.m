function [state_eci, R_rtn2eci] = rtn2eci(chief_state, state_rtn)
% rtn2eci Transform deputy satellite state from the chief's RTN frame back to ECI.
%
%   [state_eci, R_rtn2eci] = rtn2eci(chief_state, state_rtn) converts the deputy's
%   relative RTN state (position & velocity) into an absolute ECI state, given the
%   chief's absolute ECI state.
%
%   Inputs:
%     chief_state - 6×1 vector [r_c; v_c] in ECI (e.g., km, km/s)
%     state_rtn   - 6×1 vector [r_rtn; v_rtn] in the chief's RTN frame
%
%   Outputs:
%     state_eci   - 6×1 vector [r_d; v_d] in ECI
%     R_rtn2eci   - 3×3 rotation matrix from RTN → ECI
%
%   Notes:
%     Let C = [R'; T'; N'] be the matrix from ECI → RTN. Then R_rtn2eci = C'.

    % unpack chief
    rC = chief_state(1:3);
    vC = chief_state(4:6);

    % unpack deputy RTN
    r_rtn = state_rtn(1:3);
    v_rtn = state_rtn(4:6);

    % build unit vectors for chief RTN frame
    R_vec = rC / norm(rC);
    N_vec = cross(rC, vC);
    N_vec = N_vec / norm(N_vec);
    T_vec = cross(N_vec, R_vec);

    % rotation ECI → RTN
    C = [R_vec'; T_vec'; N_vec'];
    % rotation RTN → ECI
    R_rtn2eci = C';

    % recover inertial angular velocity of chief
    h = cross(rC, vC);
    omega_inertial = h / (norm(rC)^2);
    % express it in RTN
    omega_rtn = C * omega_inertial;

    % invert the position transform
    dr_eci = R_rtn2eci * r_rtn;
    rD = rC + dr_eci;

    % invert the velocity transform:
    %  v_rtn = C*(vD−vC) − ω_rtn×r_rtn
    % ⇒ C*(vD−vC) = v_rtn + ω_rtn×r_rtn
    dv = R_rtn2eci * (v_rtn + cross(omega_rtn, r_rtn));
    vD = vC + dv;

    state_eci = [rD; vD];
end
