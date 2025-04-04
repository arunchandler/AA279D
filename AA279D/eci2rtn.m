function [rv_rtn, R_rtn2eci] = eci2rtn(rv_eci)
% Convert from ECI frame to RTN frame
% Inputs:
%   rv_eci : 6x1 states vector in ECI frame [m] and [m/s]
% Outputs:
%   r_rtn : 3x1 position vector in RTN frame [m]
%   v_rtn : 3x1 velocity vector in RTN frame [m/s]
%   R_rtn2eci : 3x3 rotation matrix from RTN to ECI

r_eci = rv_eci(1:3);
v_eci = rv_eci(4:6);

R_hat = r_eci / norm(r_eci);

h = cross(r_eci, v_eci);
N_hat = h / norm(h);

T_hat = cross(N_hat, R_hat);

% Assemble the RTN-to-ECI rotation matrix
R_rtn2eci = [R_hat, T_hat, N_hat];

% Invert to get ECI-to-RTN transformation
R_eci2rtn = R_rtn2eci';

% Transform the vectors
r_rtn = R_eci2rtn * r_eci;
v_rtn = R_eci2rtn * v_eci;

rv_rtn = [r_rtn; v_rtn];

end