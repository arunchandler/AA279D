function [rv_RTN] = propagateLinearEcc( ...
    roe_qns, chief_oe, t_grid)
%PROPAGATELINEARECC  Full linear mapping for eccentric orbits
%
% Inputs:
%   roe_qns   = [δa; δλ; δe_x; δe_y; δi_x; δi_y]   (6×1)
%   chief_oe  = [a; e; i; Ω; ω; M0]              (6×1)
%   t_grid    = [t0…tN]′ times [s]              (N×1)
%   tol       = tolerance for mean2ecc solver   (scalar)
%
% Outputs:
%   r_RTN     = N×3 RTN relative positions [m]
%   v_RTN     = N×3 RTN relative velocities [m/s]

    global mu
    % unpack chief orbit
    a     = chief_oe(1);
    e     = chief_oe(2);
    i_c   = chief_oe(3);
    Omega = chief_oe(4);
    omega = chief_oe(5);
    M0    = chief_oe(6);
    nu0 = mean2true(M0, e, 1e-10);
    % constants
    n   = sqrt(mu/a^3);
    eta = sqrt(1 - e^2);

    N = numel(t_grid);
    rv_RTN = zeros(N, 6);

    for k=1:N
      t = t_grid(k);

      % 1) chief anomaly @ t
      M = M0 + n*t;
      E = mean2ecc(M, e, 1e-10);    
      f = ecc2true(E, e);         
      u = omega + f;              % true argument of latitude

      % 2) helper function
      ex   = e*cos(omega);
      ey   = e*sin(omega);
      kf   = 1 + ex*cos(u)+ ey*sin(u);
      kfp  = -ex*sin(u) + ey*cos(u);
      coti = cot(i_c);
      coss = cos(u);
      sins = sin(u);

      % 3) build B_pos (3×6)
      bxp1 = 1/kf + 1.5*kfp*(n/eta^3)*t;
      bxp2 = - kfp/eta^3;
      bxp3 =   ( ex*(kf-1)/(1+eta) - coss ) / eta^3;
      bxp4 =   ( ey*(kf-1)/(1+eta) - sins ) / eta^3;
      bxp6 =   ( kfp/eta^3 ) * coti;

      byp1 = -1.5 * kf * (n/eta^3)*t;
      byp2 =  kf/eta^3;
      byp3 = ( (1+1/kf)*sins + ey/kf + (kf/eta)*(ey/(1+eta)) ) / eta^2;
      byp4 = -( (1+1/kf)*coss + ex/kf + (kf/eta)*(ex/(1+eta)) ) / eta^2;
      byp6 = ( 1/kf - kf/eta^3 ) * coti;

      bzp5 = (1/kf)*sins;
      bzp6 =-(1/kf)*coss;

      bxp1d = kfp/2 + 1.5*kf^2*(1-kf)*(n/eta^3)*t;
      bxp2d = kf^2*(kf-1)/eta^3;
      bxp3d = kf^2/eta^3 * ( eta*sins + ey*(kf-1)/(1+eta) );
      bxp4d = -kf^2/eta^3 * ( eta*coss + ex*(kf-1)/(1+eta) );
      bxp6d = -kf^2/eta^3 * (kf-1)*coti;

      byp1d = -1.5*kf*( 1 + kf*kfp*(n/eta^3)*t );
      byp2d =  kf^2 * kfp/eta^3;
      byp3d = (1 + kf^2/eta^3)*coss + (ex*(kf/eta^2))*(1 + (kf/eta)*(1-kf)/(1+eta));
      byp4d = (1 + kf^2/eta^3)*sins + (ey*(kf/eta^2))*(1 + (kf/eta)*(1-kf)/(1+eta));
      byp6d = -(1 + kf^2/eta^3)*kfp*coti;

      bzp5d =   coss + ex;
      bzp6d = sins + ey;

      B = [bxp1, bxp2, bxp3, bxp4,   0, bxp6;
           byp1, byp2, byp3, byp4,   0, byp6;
           0,    0,    0,    0, bzp5, bzp6; 
           bxp1d, bxp2d, bxp3d, bxp4d,    0, bxp6d;
           byp1d, byp2d, byp3d, byp4d,    0, byp6d;
           0,      0,      0,      0, bzp5d, bzp6d ];

      A = [a * (eta^2)*eye(3), zeros(3); zeros(3), a * (n / eta) * eye(3)];
      tot = A * B;

      rv_RTN(k, :) = (tot * roe_qns)';
    end
end