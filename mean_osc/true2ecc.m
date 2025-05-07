% TRUE2ECC computes eccentric anomaly given true anomaly and eccentricity.
%
%   Inputs:
%     f - mean anomaly [rad]
%     e - eccentricity [-]
%
%   Outputs:
%     E - eccentric anomaly [rad]

function E = true2ecc(f, e)



f = wrapTo2Pi(f);

E = 2*atan(sqrt((1-e)/(1+e)) * tan(f/2));

% Make sure E sits in the correct semi-plane
E = wrapTo2Pi(E);

end
