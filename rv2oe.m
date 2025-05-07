function oe = rv2oe(X, mu)
    % Convert state vector to Keplerian elements [a,e,i,RAAN,omega,nu]

    r = X(1:3);
    v = X(4:6);
    mu = mu;

    h = cross(r,v);
    h_norm = norm(h);
    i = acos(h(3)/h_norm);

    K = [0;0;1];
    N = cross(K,h);
    N_norm = norm(N);
    if N_norm~=0
        RAAN = acos(N(1)/N_norm);
        if N(2)<0, RAAN = 2*pi-RAAN; end
    else
        RAAN = 0;
    end

    r_norm = norm(r);
    v_norm = norm(v);
    e_vec = (1/mu)*((v_norm^2 - mu/r_norm)*r - dot(r,v)*v);
    e     = norm(e_vec);

    if N_norm~=0
        omega = acos(dot(N,e_vec)/(N_norm*e));
        if e_vec(3)<0, omega = 2*pi-omega; end
    else
        omega = 0;
    end

    nu = acos(dot(e_vec,r)/(e*r_norm));
    if dot(r,v)<0, nu = 2*pi - nu; end

    % semi-major axis
    energy = v_norm^2/2 - mu/r_norm;
    if energy~=0
        a = -mu/(2*energy);
    else
        a = Inf;
    end

    % only return the six elements
    oe = [ a, ...
           e, ...
           wrapTo2Pi(i), ...
           wrapTo2Pi(RAAN), ...
           wrapTo2Pi(omega), ...
           wrapTo2Pi(nu) ];
end
