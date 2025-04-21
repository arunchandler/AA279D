function Phi = buildYAphi(a,e,f,f0,tau)
  global mu;
  %--- dimless functions ---
  k    = 1 + e*cos(f);
  kp   = -e*sin(f);
  dsin = cos(f) + e*cos(2*f);
  dcos = - ( sin(f) + e*sin(2*f) );
  
  Phi = zeros(6,6);
  % row 1
  Phi(1,1) = 1 + 1.5*k*kp*tau;
  Phi(1,2) = k*sin(f);
  Phi(1,3) = k*cos(f);
  % row 2
  Phi(2,1) = -1.5*k^2*tau;
  Phi(2,2) = (1 + k)*cos(f);
  Phi(2,3) = -(1 + k)*sin(f);
  Phi(2,4) = 1;
  % row 3
  Phi(3,5) = sin(f);
  Phi(3,6) = cos(f);
  % row 4
  Phi(4,1) = 1.5*(kp/k) - 1.5*e*dsin*tau;
  Phi(4,2) = dsin;
  Phi(4,3) = dcos;
  % row 5
  Phi(5,1) = -1.5*(1 + 2*k*kp*tau);
  Phi(5,2) = -2*k*sin(f);
  Phi(5,3) = e - 2*k*cos(f);
  % row 6
  Phi(6,5) = cos(f);
  Phi(6,6) = -sin(f);
end
