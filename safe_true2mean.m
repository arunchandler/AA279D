function M = safe_true2mean(f, e)
%SAFE_TRUE2MEAN  Safely compute mean anomaly given true anomaly and eccentricity
%
%   M = safe_true2mean(f, e) returns the mean anomaly, with safeguards
%   against complex or invalid values.
%
%   Inputs:
%     f - true anomaly [rad]
%     e - eccentricity [-]
%
%   Outputs:
%     M - mean anomaly [rad]

    % Check for invalid inputs
    if ~isreal(f) || ~isreal(e) || isnan(f) || isnan(e)
        fprintf('WARNING: Invalid inputs to safe_true2mean: f=%.6f, e=%.6f\n', f, e);
        M = 0;
        return;
    end
    
    % Clamp eccentricity to prevent issues
    if e >= 1.0
        fprintf('WARNING: Eccentricity >= 1.0 in safe_true2mean. Clamping to 0.9999\n');
        e = 0.9999;
    end
    
    % Use the original function
    M = true2mean(f, e);
    
    % Check output
    if ~isreal(M) || isnan(M)
        fprintf('WARNING: true2mean returned invalid result. Using 0.\n');
        M = 0;
    end
end 