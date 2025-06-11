function burns = plan_reconfiguration_dbg(chief, dep, roe_target_m, t_cur, a_chief)
% Gauss inverse with built-in debug prints  (ONLY for test harness)
% -----------------------------------------------------------------
    global Re J2 mu tol

    burns = [];

    %% ---- state conversion & requirements ------------------------
    roe_cur  = compute_roes(chief, dep)';             % dimensionless
    roe_tgt  = roe_target_m(:)' ./ a_chief;
    de_req   = roe_tgt(3:4) - roe_cur(3:4);
    di_req   = roe_tgt(5:6) - roe_cur(5:6);
    dl_req   = roe_tgt(2)   - roe_cur(2);             % raw λ gap

    %% ---- perturbed mean motion ----------------------------------
    n     = sqrt(mu/a_chief^3);
    eta   = sqrt(1 - chief(2)^2);
    P     = 3*cos(chief(3))^2 - 1;
    Q     = 5*cos(chief(3))^2 - 1;
    kappa = 0.75*J2*Re^2*sqrt(mu)/(a_chief^(7/2)*eta^4);
    nbar  = n + kappa*(eta*P + Q);

    wrap  = @(x) wrapTo2Pi(mod(x,2*pi));
    u_c   = wrap(mean2true(chief(6), chief(2), tol) + chief(5));

    %% ---- IN-PLANE pair ------------------------------------------
    if norm(de_req) > 1e-6 || abs(dl_req) > 1e-6
        xi = atan2(de_req(1), de_req(2));
        k  = 0;
        while true
            u1 = xi + k*pi; u2 = u1 + pi;
            t1 = t_cur + mod(u1 - u_c, 2*pi)/nbar;
            t2 = t_cur + mod(u2 - u_c, 2*pi)/nbar;
            if t1>t_cur && t2>t_cur, break; end, k=k+1
        end
        Delta_t = t2 - t1;

        % J2 drift compensation -----------------------------------
        gamma   = 0.5 * J2 * (Re/a_chief)^2 / eta^4;
        du_J2   = -12 * gamma * sin(2*chief(3)) * roe_cur(5) * nbar * Delta_t;
        dl_comp = dl_req - du_J2;        % compensated λ

        %%% DEBUG print block
        fprintf('---- in-plane planning ----\n');
        fprintf('raw  λ gap      : %+9.3f m\n', dl_req * a_chief);
        fprintf('du_J2 (drift)   : %+9.3f m\n', du_J2 * a_chief);
        fprintf('compensated λ   : %+9.3f m\n\n', dl_comp * a_chief);
        %%%

        de_n   = norm(de_req);
        dvr1   = n*a_chief/2 * (-dl_comp/2 + de_n);
        dvr2   = n*a_chief/2 * (-dl_comp/2 - de_n);
        max_r  = 0.5; dvr1 = sign(dvr1)*min(abs(dvr1),max_r);
                      dvr2 = sign(dvr2)*min(abs(dvr2),max_r);

        burns = [burns;
                 struct('t', t1, 'dv_rtn', [dvr1;0;0], 'tag','R1', 'type','reconfig');
                 struct('t', t2, 'dv_rtn', [dvr2;0;0], 'tag','R2', 'type','reconfig')];
    end

    %% ---- OUT-OF-PLANE single pulse ------------------------------
    if norm(di_req) > 1e-6
        if isempty(burns)
            t_ref = t_cur;                  % no R1/R2 in the queue
        else
            t_ref = max([t_cur; [burns.t]']);
        end
        k=0; 
        while true
            uN = atan2(di_req(2), di_req(1)) + k*pi;
            tn = t_ref + mod(uN - u_c, 2*pi) / nbar;
            if tn>t_ref, break; end, k=k+1;
        end
        dvn  = n*a_chief*norm(di_req);
        dvn  = sign(dvn)*min(abs(dvn),0.5);

        burns = [burns;
                 struct('t', tn, 'dv_rtn', [0;0;dvn],'tag','N','type','reconfig')];
    end
end
