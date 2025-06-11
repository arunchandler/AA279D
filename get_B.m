function B = get_B(oe)
%GET_B Compute the Gauss variational equation B matrix for a chief orbit
%   B = GET_B(oe) returns the 6×3 B matrix mapping perturbing accelerations
%   [u_r; u_t; u_n] to orbital element rates. Input:
%   oe = [a; e; i; omega; f] (semi-major axis, eccentricity, inclination,
%        argument of perigee, true anomaly).
%   Uses global mu (gravitational parameter).

global mu;

% Unpack orbital elements
a_c     = oe(1);
e_c     = oe(2);
i_c     = oe(3);
omega_c = oe(4);
f_c     = oe(5);

% Auxiliary parameters
eta_c = sqrt(1 - e_c^2);
n_c   = sqrt(mu / a_c^3);
ex    = e_c * cos(omega_c);
ey    = e_c * sin(omega_c);
denom = 1 + e_c * cos(f_c);

% Build unscaled B-tilde matrix
Btilde = zeros(6,3);

Btilde(1,1) = 2*e_c*sin(f_c)/eta_c;
Btilde(1,2) = 2*(1 + e_c*cos(f_c))/eta_c;

Btilde(2,1) = -2*eta_c^2/denom;

Btilde(3,1) = eta_c*sin(omega_c + f_c);
Btilde(3,2) = eta_c*((2 + e_c*cos(f_c))*cos(omega_c + f_c) + ex)/denom;
Btilde(3,3) = eta_c*ey*sin(omega_c + f_c)/(tan(i_c)*denom);

Btilde(4,1) = -eta_c*cos(omega_c + f_c);
Btilde(4,2) = eta_c*((2 + e_c*cos(f_c))*sin(omega_c + f_c) + ey)/denom;
Btilde(4,3) = -eta_c*ex*sin(omega_c + f_c)/(tan(i_c)*denom);

Btilde(5,3) = eta_c*cos(omega_c + f_c)/denom;
Btilde(6,3) = eta_c*sin(omega_c + f_c)/denom;

% Scale by 1/(a * n)
B = Btilde/(a_c * n_c);
end
