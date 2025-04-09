% ECC2TRUE computes true anomaly given eccentric anomaly and eccentricity.
%
%   Inputs:
%     E - eccentric anomaly [rad]
%     e - eccentricity [-]
%
%   Outputs:
%     f - true anomaly [rad]

function f = ecc2true(E, e)

E = wrapTo2Pi(E);

f = atan2(sin(E)*sqrt(1-e^2),cos(E)-e);

end
