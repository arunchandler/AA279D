function [deputyOE, success] = safe_qns2oe(chiefOE, dqns, max_ecc)
%SAFE_QNS2OE  Safely convert chief elements + QNS-relative elements to deputy OE
%
%  [deputyOE, success] = safe_qns2oe(chiefOE, dqns, max_ecc) returns the deputy's 
%  classical orbital elements [a; e; i; Omega; omega; M] given
%    chiefOE = [a_c; e_c; i_c; Omega_c; omega_c; M_c]
%    dqns    = a_c*[delta_a; delta_lambda; delta_ex; delta_ey; delta_ix; delta_iy]
%    max_ecc = maximum allowed eccentricity (default: 0.9999)
%
%  If the computed eccentricity exceeds max_ecc, it clamps the eccentricity
%  and returns success = false.

    if nargin < 3
        max_ecc = 0.9999;
    end

    % unpack
    a_c     = chiefOE(1);
    e_c     = chiefOE(2);
    i_c     = chiefOE(3);
    Omega_c = chiefOE(4);
    omega_c = chiefOE(5);
    M_c     = chiefOE(6);

    da      = dqns(1)/a_c;
    dlambda = dqns(2)/a_c;
    dex     = dqns(3)/a_c;
    dey     = dqns(4)/a_c;
    dix     = dqns(5)/a_c;
    diy     = dqns(6)/a_c;

    % semi‐major axis
    a_d = a_c*(1 + da);

    % inclination
    i_d = i_c + dix;

    % RAAN
    %   di_y = (RAAN_d - RAAN_c) * sin(i_c)  →  RAAN_d = RAAN_c + di_y/sin(i_c)
    RAAN_d = Omega_c + diy/sin(i_c);

    % eccentricity vector in the orbital plane
    evec_c = [e_c*cos(omega_c); e_c*sin(omega_c)];
    evec_d = evec_c + [dex; dey];
    e_d    = hypot(evec_d(1), evec_d(2));           % = sqrt(x^2+y^2)
    
    % Check if eccentricity needs clamping
    if e_d >= max_ecc
        fprintf('WARNING: Computed eccentricity = %.6f >= %.4f. Clamping to %.4f\n', ...
                e_d, max_ecc, max_ecc);
        % Scale the eccentricity vector to maintain direction but limit magnitude
        evec_d = evec_d * (max_ecc / e_d);
        e_d = max_ecc;
        success = false;
    else
        success = true;
    end
    
    omega_d = atan2(evec_d(2), evec_d(1));           % true anomaly argument

    % mean longitude difference gives M_d
    %   dlambda = (M_d+omega_d)-(M_c+omega_c)+(RAAN_d-RAAN_c)*cos(i_c)
    M_plus_omega = dlambda + (M_c + omega_c) - (RAAN_d - Omega_c)*cos(i_c);
    M_d          = M_plus_omega - omega_d;

    % wrap angles into [0,2π)
    RAAN_d = mod(RAAN_d,2*pi);
    omega_d = mod(omega_d,2*pi);
    M_d = mod(M_d,2*pi);

    deputyOE = [a_d; e_d; i_d; RAAN_d; omega_d; M_d];
end 