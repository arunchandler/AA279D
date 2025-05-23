function params = rv2oe(X, mu)
    % Convert state vector (position and velocity) to Keplerian orbital elements
    % X: State vector [rx; ry; rz; vx; vy; vz]
    % mu: Gravitational parameter of the central body
    % params: Array [a, e, i, RAAN, omega, nu, energy, h, e_vec] where h is
    % 1x3 and e_vec is 1x2

    r = X(1:3);
    v = X(4:6);

    r_norm = norm(r);
    v_norm = norm(v);

    h = cross(r, v);
    h_norm = norm(h);

    % Inclination
    i = wrapTo2Pi(acos(h(3) / h_norm));

    % Node vector
    K = [0; 0; 1];
    N = cross(K, h);
    N_norm = norm(N);

    % Right Ascension of the Ascending Node (RAAN)
    if N_norm ~= 0
        RAAN = wrapTo2Pi(acos(N(1) / N_norm));
        if N(2) < 0
            RAAN = wrapTo2Pi(2*pi - RAAN);
        end
    else
        RAAN = 0;
    end

    % Eccentricity vector
    e_vec = (1/mu) * ( (v_norm^2 - mu/r_norm) * r - dot(r, v) * v );
    e = norm(e_vec);

    % Argument of Periapsis (omega)
    if N_norm ~= 0
        omega = wrapTo2Pi(acos(dot(N, e_vec) / (N_norm * e)));
        if e_vec(3) < 0
            omega = wrapTo2Pi(2*pi - omega);
        end
    else
        omega = 0;
    end

    % True Anomaly (nu)
    nu = wrapTo2Pi(acos(dot(e_vec, r) / (e * r_norm)));
    if dot(r, v) < 0
        nu = wrapTo2Pi(2*pi - nu);
    end

    % Semi-major axis (a)
    energy = (v_norm^2 / 2) - (mu / r_norm);
    if energy ~= 0
        a = -mu / (2 * energy);
    else
        a = Inf;
    end

    e_vec_perifocal = [e*cos(omega), e*sin(omega)];

    params = [a, e, wrapTo2Pi(i), wrapTo2Pi(RAAN), wrapTo2Pi(omega), wrapTo2Pi(nu), energy, h(:).', e_vec_perifocal];
end
