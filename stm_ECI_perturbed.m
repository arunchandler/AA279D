function phi = stm_ECI_perturbed(tau,rv)
%Calculating STM for state propagation
  % build Jacobian F at the current "reference" rv
  r      = rv(1:3);   x=r(1); y=r(2); z=r(3);
  r_norm = norm(r);   r2 = r_norm^2;   r5=r_norm^5;  r7=r_norm^7;
  K      = 1.5 * J2 * mu * Re^2;

  % ∂a/∂r (exactly as in your linearized function)
  daxdx = -mu*(1/r_norm^3 - 3*x^2/r5) ...
          - K*(1/r5*(1-5*z^2/r2) - 5*x^2/r7*(3-5*z^2/r2));
  daxdy =  3*mu*x*y/r5 ...
          - K*(-5*x*y/r7*(3-5*z^2/r2));
  daxdz =  3*mu*x*z/r5 ...
          - K*(-10*x*z/r7*(1 - z^2/r2));

  daydx =  3*mu*y*x/r5 ...
          - K*(-5*y*x/r7*(3-5*z^2/r2));
  daydy = -mu*(1/r_norm^3 - 3*y^2/r5) ...
          - K*(1/r5*(1-5*z^2/r2) - 5*y^2/r7*(3-5*z^2/r2));
  daydz =  3*mu*y*z/r5 ...
          - K*(-10*y*z/r7*(1 - z^2/r2));

  dazdx =  3*mu*z*x/r5 ...
          - K*(-5*z*x/r7*(3-5*z^2/r2));
  dazdy =  3*mu*z*y/r5 ...
          - K*(-5*z*y/r7*(3-5*z^2/r2));
  dazdz = -mu*(1/r_norm^3 - 3*z^2/r5) ...
          - K*(3/r5 - 5*z^2/r7*(6-5*z^2/r2));

  F = zeros(6);
  F(1:3,4:6)    = eye(3);
  F(4,1:3)      = [daxdx, daxdy, daxdz];
  F(5,1:3)      = [daydx, daydy, daydz];
  F(6,1:3)      = [dazdx, dazdy, dazdz];
end