function params = rv2oe(X, mu)
    % RV2OE  Convert state vector to classical orbital elements robustly
    %   params = rv2oe(X, mu) returns
    %     [a, e, i, RAAN, omega, nu, energy, h_x, h_y, h_z, e_pf_x, e_pf_y]
    %
    %  — clamps acos‐arguments to [-1,1]
    %  — forces 1<e<1+tol → e=1 for near‐parabolic orbits

    % unpack
    r = X(1:3);  v = X(4:6);
    r_norm = norm(r);
    v_norm = norm(v);

    % angular momentum
    h_vec  = cross(r,v);
    h_norm = norm(h_vec);

    % inclination
    i = acos(clamp(h_vec(3)/h_norm, -1, 1));

    % node line
    K     = [0;0;1];
    N_vec = cross(K, h_vec);
    N_norm = norm(N_vec);

    % RAAN
    if N_norm > eps
        RAAN = acos(clamp(N_vec(1)/N_norm, -1, 1));
        if N_vec(2) < 0
            RAAN = 2*pi - RAAN;
        end
    else
        RAAN = 0;
    end

    % eccentricity vector & magnitude
    e_vec = (1/mu)*((v_norm^2 - mu/r_norm)*r - dot(r,v)*v);
    e = norm(e_vec);

    % if it’s supposed to be elliptical but numerics pushed e>1
    tol_e = 1e-10;
    if e>1 && e<1+tol_e
        e = 1;
        e_vec = e_vec / norm(e_vec);
    end

    % argument of periapsis
    if N_norm>eps && e>tol_e
        omega = acos(clamp(dot(N_vec,e_vec)/(N_norm*e), -1, 1));
        if e_vec(3) < 0
            omega = 2*pi - omega;
        end
    else
        omega = 0;
    end

    % true anomaly
    if e>tol_e
        nu = acos(clamp(dot(e_vec,r)/(e*r_norm), -1, 1));
        if dot(r,v) < 0
            nu = 2*pi - nu;
        end
    else
        nu = 0;
    end

    % specific energy & semi‐major axis
    energy = v_norm^2/2 - mu/r_norm;
    if abs(energy)>eps
        a = -mu/(2*energy);
    else
        a = Inf;
    end

    % eject tiny negative zeros
    i    = wrapTo2Pi(i);
    RAAN = wrapTo2Pi(RAAN);
    omega= wrapTo2Pi(omega);
    nu   = wrapTo2Pi(nu);

    % build output
    h_pf = h_vec.';                              % [h_x h_y h_z]
    e_pf = [e*cos(omega), e*sin(omega)];         % eccentricity in perifocal

    params = [a, e, i, RAAN, omega, nu, energy, h_pf, e_pf];
end

function y = clamp(x, lo, hi)
    % simple helper to bound x in [lo, hi]
    y = min( max(x, lo), hi );
end
