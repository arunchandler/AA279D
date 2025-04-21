function Phi = buildHCWphi(a,t)
    global mu;
    n = sqrt(mu/a^3);
    nt = n*t;
    an = a*n;
    
    mat1 = [a, 0, 0, 0, 0, 0;
            0, a, 0, 0, 0, 0;
            0, 0, a, 0, 0, 0;
            0, 0, 0, an, 0, 0;
            0, 0, 0, 0, an, 0;
            0, 0, 0, 0, 0, an];
    
    mat2 = [1,       sin(nt),    cos(nt),    0, 0,       0;
            -3/2*nt, 2*cos(nt),  -2*sin(nt), 1, 0,       0;
            0,       0,          0,          0, sin(nt), cos(nt);
            0,       cos(nt),    -sin(nt),   0, 0,       0;
            -3/2,    -2*sin(nt), -2*cos(nt), 0, 0,       0;
            0,       0,          0,          0, cos(nt), -sin(nt)];
    
    Phi = mat1 * mat2;

end