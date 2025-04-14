function [state_rtn, R_eci2rtn] = eci2rtn(chief_state, deputy_state)
% eci2rtn Transform deputy satellite state from ECI to the chief's RTN frame.
%
%   state_rtn = eci2rtn(chief_state, deputy_state) converts the deputy's state 
%   (position and velocity) from the Earth-Centered Inertial (ECI) frame to the 
%   Radial-Transverse-Normal (RTN) frame of the chief satellite.
%
%   Inputs:
%       chief_state  - 6x1 vector [r_c; v_c] where:
%                        r_c: 3x1 position vector in ECI (e.g., in km)
%                        v_c: 3x1 velocity vector in ECI (e.g., in km/s)
%       deputy_state - 6x1 vector [r_d; v_d] with similar structure for the deputy.
%
%   Output:
%       state_rtn - 6x1 vector [r_rtn; v_rtn] where:
%                     r_rtn: 3x1 relative position of the deputy expressed in the 
%                            chief's RTN frame.
%                     v_rtn: 3x1 relative velocity of the deputy in the RTN frame,
%                            which includes the correction for the rotating frame.
%
%   The RTN frame is defined as:
%       R = r_c / ||r_c|| (Radial direction)
%       N = (r_c x v_c) / ||r_c x v_c|| (Normal to the orbit)
%       T = N x R (Transverse direction)

    % Chief
    rC = chief_state(1:3);
    vC = chief_state(4:6);
    
    % Deputy
    rD = deputy_state(1:3);
    vD = deputy_state(4:6);
    
    %chief's RTN coordinate system
    R_vec = rC / norm(rC);
    N_vec = cross(rC, vC);
    N_vec = N_vec / norm(N_vec);
    T_vec = cross(N_vec, R_vec);
    
    %Rotation matrix
    C = [R_vec'; T_vec'; N_vec'];
    
    dr = rD - rC;
    dv = vD - vC;
    
    r_rtn = C * dr;
    
    % Compute the angular velocity of the chief's orbit
    h = cross(rC, vC);
    omega_inertial = h / (norm(rC)^2);
    
    % Express the angular velocity in the RTN frame.
    omega_rtn = C * omega_inertial;
    
    % Transform the relative velocity to the RTN frame
    v_rtn = C * dv - cross(omega_rtn, r_rtn);
    
    state_rtn = [r_rtn; v_rtn];
    R_eci2rtn = C;
end
