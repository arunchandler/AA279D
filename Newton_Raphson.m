function E = Newton_Raphson(M,e,epsilon) % takes radians, outputs radians
    E = M;

    while true
        fE = E - e*sin(E) - M;
        fprimeE = 1 - e*cos(E);
        delta = -fE/fprimeE;
        E_next = E+delta;
        if abs(delta)<epsilon
            break;
        end
    
        E = E_next;
    end
end