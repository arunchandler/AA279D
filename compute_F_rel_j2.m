function F_rel = compute_F_rel_j2(delta_r, r_chief)
  % delta_r: 3×1 relative position (r_dep - r_chief)
  % r_chief: 3×1 chief absolute ECI position
  % returns 6×6 Jacobian for [delta_r; delta_v]

  % absolute deputy position
  r_dep = r_chief + delta_r;

  % compute the 3×3 accel‐Jacobian at r_dep
  A_dep = compute_acc_jacobian_j2(r_dep);

  F_rel = zeros(6);
  F_rel(1:3,4:6) = eye(3);
  F_rel(4:6,1:3) = A_dep;
end

function A = compute_acc_jacobian_j2(r)
  global mu J2 Re
  x = r(1); y = r(2); z = r(3);
  R = norm(r); R2=R^2; R5=R^5; R7=R^7;
  K = 1.5*J2*mu*Re^2;
  z2 = z^2;

  % fill A(1,1), A(1,2), A(1,3) exactly as in your compute_F_j2
  A = zeros(3);
  A(1,1) = -mu*(1/R^3 - 3*x^2/R5) ...
           -K*((1/R5)*(1-5*z2/R2) - 5*x^2/R7*(3-5*z2/R2));
  A(1,2) =  3*mu*x*y/R5 ...
           -K*(-5*x*y/R7*(3-5*z2/R2));
  A(1,3) =  3*mu*x*z/R5 ...
           -K*(-10*x*z/R7*(1 - z2/R2));

  % and similarly A(2,1)=A(1,2), A(2,2), A(2,3), A(3,1)=A(1,3),
  % A(3,2)=A(2,3), A(3,3) as above…
  A(2,1) = A(1,2);
  A(2,2) = -mu*(1/R^3 - 3*y^2/R5) ...
           -K*((1/R5)*(1-5*z2/R2) - 5*y^2/R7*(3-5*z2/R2));
  A(2,3) =  3*mu*y*z/R5 ...
           -K*(-10*y*z/R7*(1 - z2/R2));
  A(3,1) = A(1,3);
  A(3,2) = A(2,3);
  A(3,3) = -mu*(1/R^3 - 3*z2/R5) ...
           -K*((3/R5) - 5*z2/R7*(6-5*z2/R2));
end