% oe = [a; e; i; RAAN; ω; M]
% qnsoe = [a; ex; ey; i; RAAN; u]
function qnsoe = classical2qns(oe)
    % Unpack classical elements
    a     = oe(1);    % semi‑major axis
    e     = oe(2);    % eccentricity
    i     = oe(3);    % inclination
    RAAN  = oe(4);    % right ascension of ascending node
    omega = oe(5);    % argument of perigee
    M     = oe(6);    % mean anomaly

    % Compute quasi‑nonsingular components
    ex = e .* cos(omega);
    ey = e .* sin(omega);
    u  = M + omega;           % mean argument of latitude

    % (Optional) wrap angles into [0, 2π)
    RAAN = mod(RAAN, 2*pi);
    u    = mod(u,    2*pi);

    % Assemble output
    qnsoe = [a; ex; ey; i; RAAN; u];
end
