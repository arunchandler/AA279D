function [oe, success] = safe_rv2oe(rv, mu, max_ecc)
%SAFE_RV2OE  Safely convert position-velocity to orbital elements
%
%   [oe, success] = safe_rv2oe(rv, mu, max_ecc) returns the orbital elements
%   and a success flag. If eccentricity exceeds max_ecc, it clamps the
%   eccentricity and returns success = false.
%
%   Inputs:
%     rv     - [r; v] state vector (6x1)
%     mu     - gravitational parameter
%     max_ecc - maximum allowed eccentricity (default: 0.9999)
%
%   Outputs:
%     oe     - [a, e, i, RAAN, omega, nu, energy, h_pf(3), e_pf(2)]
%     success - true if no eccentricity clamping was needed

    if nargin < 3
        max_ecc = 0.9999;
    end
    
    % Call the original rv2oe function
    oe = rv2oe(rv, mu);
    
    % Check if eccentricity needs clamping
    if oe(2) >= max_ecc
        fprintf('WARNING: Eccentricity = %.6f >= %.4f. Clamping to %.4f\n', ...
                oe(2), max_ecc, max_ecc);
        oe(2) = max_ecc;
        success = false;
    else
        success = true;
    end
end 