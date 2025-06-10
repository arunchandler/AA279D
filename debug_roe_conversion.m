function debug_roe_conversion(chief_oe, deputy_oe, rel_qns_init)
%DEBUG_ROE_CONVERSION  Debug ROE conversion process
%
%   debug_roe_conversion(chief_oe, deputy_oe, rel_qns_init) analyzes the
%   ROE conversion process to identify potential issues.

    fprintf('\n=== ROE Conversion Debug ===\n');
    
    % Unpack chief elements
    a_c = chief_oe(1);
    e_c = chief_oe(2);
    i_c = chief_oe(3);
    RAAN_c = chief_oe(4);
    omega_c = chief_oe(5);
    M_c = chief_oe(6);
    
    % Unpack deputy elements
    a_d = deputy_oe(1);
    e_d = deputy_oe(2);
    i_d = deputy_oe(3);
    RAAN_d = deputy_oe(4);
    omega_d = deputy_oe(5);
    M_d = deputy_oe(6);
    
    fprintf('Chief OE:  a=%.2f, e=%.6f, i=%.4f°, RAAN=%.4f°, ω=%.4f°, M=%.4f°\n', ...
            a_c, e_c, rad2deg(i_c), rad2deg(RAAN_c), rad2deg(omega_c), rad2deg(M_c));
    fprintf('Deputy OE: a=%.2f, e=%.6f, i=%.4f°, RAAN=%.4f°, ω=%.4f°, M=%.4f°\n', ...
            a_d, e_d, rad2deg(i_d), rad2deg(RAAN_d), rad2deg(omega_d), rad2deg(M_d));
    
    % Compute ROEs using compute_roes (dimensionless)
    roe_hat = compute_roes(chief_oe, deputy_oe)';
    fprintf('\nROEs (dimensionless):\n');
    fprintf('  δa = %.6f\n', roe_hat(1));
    fprintf('  δλ = %.6f\n', roe_hat(2));
    fprintf('  δex = %.6f\n', roe_hat(3));
    fprintf('  δey = %.6f\n', roe_hat(4));
    fprintf('  δix = %.6f\n', roe_hat(5));
    fprintf('  δiy = %.6f\n', roe_hat(6));
    
    % Convert to meters
    roe_m = a_c * roe_hat;
    fprintf('\nROEs (meters):\n');
    fprintf('  a·δa = %.2f m\n', roe_m(1));
    fprintf('  a·δλ = %.2f m\n', roe_m(2));
    fprintf('  a·δex = %.2f m\n', roe_m(3));
    fprintf('  a·δey = %.2f m\n', roe_m(4));
    fprintf('  a·δix = %.2f m\n', roe_m(5));
    fprintf('  a·δiy = %.2f m\n', roe_m(6));
    
    % Compare with input rel_qns_init
    if nargin >= 3
        fprintf('\nInput rel_qns_init (meters):\n');
        fprintf('  [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]\n', rel_qns_init);
        
        fprintf('\nDifference (computed - input):\n');
        diff_m = roe_m - rel_qns_init;
        fprintf('  [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f]\n', diff_m);
    end
    
    % Check for potential issues
    fprintf('\nPotential Issues:\n');
    
    % Check eccentricity
    if e_d >= 1.0
        fprintf('  WARNING: Deputy eccentricity >= 1.0 (e = %.6f)\n', e_d);
    elseif e_d > 0.99
        fprintf('  WARNING: Deputy eccentricity very high (e = %.6f)\n', e_d);
    end
    
    % Check semi-major axis
    if a_d <= 0
        fprintf('  ERROR: Deputy semi-major axis <= 0 (a = %.2f)\n', a_d);
    end
    
    % Check inclination
    if i_d < 0 || i_d > pi
        fprintf('  WARNING: Deputy inclination out of range (i = %.4f°)\n', rad2deg(i_d));
    end
    
    % Check ROE magnitudes
    if norm(roe_hat(3:4)) > 0.1
        fprintf('  WARNING: Large relative eccentricity (||δe|| = %.6f)\n', norm(roe_hat(3:4)));
    end
    
    if norm(roe_hat(5:6)) > 0.1
        fprintf('  WARNING: Large relative inclination (||δi|| = %.6f)\n', norm(roe_hat(5:6)));
    end
    
    fprintf('=== End Debug ===\n\n');
end 