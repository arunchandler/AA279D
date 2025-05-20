function Phi = buildYAphi(a, e, f, f0, tau)
    % YAMANAKA-ANKERSEN transition matrix at true anomaly f
    %   Phi = A·B exactly as in the reference
    global mu;

    % orbit parameters
    n   = sqrt(mu/a^3);
    eta = sqrt(1 - e^2);
    
    % dimensionless B-matrix entries
    k   = 1 + e*cos(f);
    kp  = -e*sin(f);

    psi_x1 = 1/k + 1.5*kp * tau;
    psi_x2 =         sin(f);
    psi_x3 =         cos(f);

    psi_y1 = -1.5*k * tau;
    psi_y2 = (1 + 1/k)*cos(f);
    psi_y3 =-(1 + 1/k)*sin(f);
    psi_y4 = 1/k;

    psi_z5 = (1/k)*sin(f);
    psi_z6 = (1/k)*cos(f);

    psi_xd1 = (kp/2) - 1.5*(k^2)*(k-1)*tau;
    psi_xd2 =  k^2*cos(f);
    psi_xd3 = -k^2*sin(f);

    psi_yd1 = -1.5*(k + k^2*kp*tau);
    psi_yd2 = -(k^2 + 1)*sin(f);
    psi_yd3 = -e - (k^2 + 1)*cos(f);
    psi_yd4 = -kp;

    psi_zd5 =  e + cos(f);
    psi_zd6 = -sin(f);

    % assemble B
    B = zeros(6,6);
    B(1,1:3) = [psi_x1, psi_x2, psi_x3];
    B(2,1:4) = [psi_y1, psi_y2, psi_y3, psi_y4];
    B(3,5:6) = [psi_z5, psi_z6];
    B(4,1:3) = [psi_xd1, psi_xd2, psi_xd3];
    B(5,1:4) = [psi_yd1, psi_yd2, psi_yd3, psi_yd4];
    B(6,5:6) = [psi_zd5, psi_zd6];

    % scale by A = diag([a·η², a·η², a·η², a·n/η, a·n/η, a·n/η])
    A = [a*eta^2*eye(3),           zeros(3);
             zeros(3), a*(n/eta)*eye(3)];

    % final YA transition matrix
    Phi = A * B;
end