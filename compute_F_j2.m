function F = compute_F_j2(r)
  % returns the 6×6 Jacobian ∂[v; a]/∂[r; v] for two‑body + J2
  global mu J2 Re

  x = r(1); y = r(2); z = r(3);
  R = norm(r);
  R2 = R^2; R5 = R^5; R7 = R^7;
  K = 1.5 * J2 * mu * Re^2;

  % precompute
  z2 = z^2;

  % ∂a/∂r components
  daxdx = -mu*(1/R^3 - 3*x^2/R5) ...
          - K*( (1/R5)*(1-5*z2/R2) - 5*x^2/R7*(3-5*z2/R2) );
  daxdy =  3*mu*x*y/R5 ...
          - K*( -5*x*y/R7*(3-5*z2/R2) );
  daxdz =  3*mu*x*z/R5 ...
          - K*( -10*x*z/R7*(1 - z2/R2) );

  daydx = daxdy;
  daydy = -mu*(1/R^3 - 3*y^2/R5) ...
          - K*( (1/R5)*(1-5*z2/R2) - 5*y^2/R7*(3-5*z2/R2) );
  daydz =  3*mu*y*z/R5 ...
          - K*( -10*y*z/R7*(1 - z2/R2) );

  dazdx =  3*mu*z*x/R5 ...
          - K*( -5*z*x/R7*(3-5*z2/R2) );
  dazdy =  3*mu*z*y/R5 ...
          - K*( -5*z*y/R7*(3-5*z2/R2) );
  dazdz = -mu*(1/R^3 - 3*z^2/R5) ...
          - K*( (3/R5) - 5*z^2/R7*(6-5*z2/R2) );

  % assemble F = [ 0₃  I₃ ; ∂a/∂r  0₃ ]
  F = zeros(6);
  F(1:3,4:6) = eye(3);
  F(4,1:3) = [daxdx, daxdy, daxdz];
  F(5,1:3) = [daydx, daydy, daydz];
  F(6,1:3) = [dazdx, dazdy, dazdz];
end