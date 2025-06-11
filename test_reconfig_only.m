% ================================================================
%  reconfig_only.m   –  mission driver *without* formation keeping
% ================================================================
clc; clear; close all;
addpath mean_osc                 % helper folder
format long g

%% ---------- constants, helpers ----------
global tol Re J2 mu s_d
tol = 1e-10;
Re  = 6378137;                    J2 = 1.082626e-3;
mu  = 3.986004418e14;             s_d = 86400;

vrow = @(x) x(:).';               % force row
wrap = @(x) wrapTo2Pi(mod(x,2*pi));

%% ---------- chief initial orbit ----------
a_TSX    = 6886536.686;           e_TSX = 0.0001264;
i_TSX    = deg2rad(97.4453);
RAAN_TSX = deg2rad(351.0108);
omega_TSX= deg2rad(101.2452);
M_TSX    = deg2rad(11.6520);
TSX_init = [a_TSX e_TSX i_TSX RAAN_TSX omega_TSX M_TSX];

%% ---------- mission scenarios (reconfigs only) ----------
scenario = input(['Select scenario:\n' ...
                  '  1: Single reconfiguration (M-D2 → M-D3)\n' ...
                  '  2: Two reconfigurations (D2→D3→D1)\n' ...
                  'Your choice: ']);

switch scenario
    case 1      %  ▸ coast 12 h ▸ reconfig ▸ coast remainder
        phases = [struct('type','coast',    'start',0,'dur',0.5*s_d, ...
                         'qns',[0 0 0 300 0 400]);
                  struct('type','reconfig', 'start',0.5*s_d,'dur',0.1*s_d, ...
                         'qns',[0 0 0 300 0 500]);
                  struct('type','coast',    'start',0.6*s_d,'dur',0.4*s_d, ...
                         'qns',[0 0 0 300 0 500])];
        sim_days = 1;
    case 2      %  ▸ coast ▸ Rcfg1 ▸ coast ▸ Rcfg2 ▸ coast
        phases = [struct('type','coast',    'start',0,'dur',1*s_d, ...
                         'qns',[0 0 0 300 0 400]);
                  struct('type','reconfig', 'start',1*s_d,'dur',0.1*s_d, ...
                         'qns',[0 0 0 300 0 500]);
                  struct('type','coast',    'start',1.1*s_d,'dur',0.9*s_d, ...
                         'qns',[0 0 0 300 0 500]);
                  struct('type','reconfig', 'start',2*s_d,'dur',0.1*s_d, ...
                         'qns',[0 0 0 500 0 300]);
                  struct('type','coast',    'start',2.1*s_d,'dur',0.9*s_d, ...
                         'qns',[0 0 0 500 0 300])];
        sim_days = 3;
    otherwise
        error('Invalid scenario');
end

%% ---------- initialise deputy ----------
[TDX_init,~] = safe_qns2oe(TSX_init, phases(1).qns');

%% ---------- integration grid ----------
n     = sqrt(mu/a_TSX^3);
t_end = sim_days * s_d;
t_grid= linspace(0, t_end, 10001).';
dt    = t_grid(2) - t_grid(1);

%% ---------- state history arrays ----------
N  = numel(t_grid);
TSX = zeros(N,6);  TSX(1,:) = TSX_init;
TDX = zeros(N,6);  TDX(1,:) = TDX_init;
rel = zeros(N,6);  rel(1,:) = a_TSX * compute_roes(TSX_init,TDX_init)';

%% ---------- burn bookkeeping ----------
pending   = [];          % burn queue
burn_t    = [];          % executed epochs
burn_dv   = [];          % executed |Δv|
burn_vecs = [];          % executed RTN components
cum_dv    = 0;

%% ---------- MAIN LOOP ----------
phase_idx = 1;           % current phase
for k = 2:N
    t_prev = t_grid(k-1);       t_cur = t_grid(k);

    % 1) propagate chief & deputy (row vectors!)
    [~, c_tmp] = ode4(@compute_rates_GVE_J2,[t_prev t_cur]',TSX(k-1,:),dt);
    [~, d_tmp] = ode4(@compute_rates_GVE_J2,[t_prev t_cur]',TDX(k-1,:),dt);
    chief = vrow(c_tmp(end,:));           dep   = vrow(d_tmp(end,:));

    % 2) execute queued burns whose epoch ≤ t_cur
    while ~isempty(pending) && pending(1).t <= t_cur + 1e-9
        b = pending(1);

        % rewind both spacecraft to the exact burn epoch
        [~,c_tmp] = ode4(@compute_rates_GVE_J2,[t_prev b.t]',chief,dt/10);
        [~,d_tmp] = ode4(@compute_rates_GVE_J2,[t_prev b.t]',dep,  dt/10);
        chief = vrow(c_tmp(end,:));       dep   = vrow(d_tmp(end,:));

        % impulse on deputy
        rv = oe2rv(dep,mu);   dv = dv_RTN2ECI(rv(1:3),rv(4:6),b.dv_rtn);
        dep = vrow(kick(dep,dv));

        % record bookkeeping
        burn_t    = [burn_t; b.t];
        burn_dv   = [burn_dv; norm(dv)];
        burn_vecs = [burn_vecs; b.dv_rtn.'];
        cum_dv    = cum_dv + norm(dv);
        fprintf('Executed %s at t=%8.1f s |Δv|=%6.2f mm/s\n', ...
                b.tag, b.t, norm(dv)*1e3);

        % propagate to current grid point again
        [~,c_tmp] = ode4(@compute_rates_GVE_J2,[b.t t_cur]',chief,dt/10);
        [~,d_tmp] = ode4(@compute_rates_GVE_J2,[b.t t_cur]',dep,  dt/10);
        chief = vrow(c_tmp(end,:));   dep = vrow(d_tmp(end,:));

        pending(1) = [];      t_prev = b.t;   % pop burn
    end

    % 3) log states
    TSX(k,:) = chief;   TSX(k,6) = wrap(TSX(k,6));
    TDX(k,:) = dep;     TDX(k,6) = wrap(TDX(k,6));
    rel(k,:) = a_TSX * compute_roes(chief,dep)';

    % 4) advance phase pointer (simple LIN search, ok for few phases)
    if phase_idx < numel(phases) && ...
       t_cur >= phases(phase_idx).start + phases(phase_idx).dur
        phase_idx = phase_idx + 1;
    end
    phase = phases(phase_idx);

    % 5) if we just entered a re-config phase → plan burns once
    if strcmp(phase.type,'reconfig') && isempty(pending)
        burns = plan_reconfiguration(chief,dep,phase.qns',t_cur,a_TSX);
        pending = [pending; burns];
        fprintf('\n=== Reconfiguration at t=%.1f s ===\n',t_cur);
        fprintf('Target ROE: [%g %g %g %g %g %g] m\n',phase.qns);
        fprintf('Planned %d burns\n',numel(burns));
    end
end

%% ---------- summary ----------
fprintf('\n=== SUMMARY ===\n');
fprintf('Total burns   : %d\n',numel(burn_t));
fprintf('Total Δv      : %.3f m/s (%.3f m/s per day)\n', ...
        cum_dv, cum_dv/sim_days);

%% ---------- plots (optional) ----------
% … (reuse your ROE, burn-timeline, cumulative-Δv plots if desired)
%   rel(:,2) will show immediately whether λ hits the target.

% -------------------------------------------------------------------------
%  ↓↓↓  only one manoeuvre planner remains  ↓↓↓
% -------------------------------------------------------------------------
function burns = plan_reconfiguration(chief,dep,roe_tgt_m,t_cur,a)
% Gauss‐inverse, with J₂ λ-drift compensation (rows everywhere)
    global Re J2 mu tol
    burns = [];

    roe_cur = compute_roes(chief,dep)';          % dimensionless
    roe_tgt = roe_tgt_m(:)' ./ a;

    de = roe_tgt(3:4) - roe_cur(3:4);
    di = roe_tgt(5:6) - roe_cur(5:6);
    dl = roe_tgt(2)   - roe_cur(2);

    % perturbed mean motion
    n     = sqrt(mu/a^3);
    eta   = sqrt(1-chief(2)^2);
    P     = 3*cos(chief(3))^2 - 1;
    Q     = 5*cos(chief(3))^2 - 1;
    kappa = 0.75*J2*Re^2*sqrt(mu)/(a^(7/2)*eta^4);
    nbar  = n + kappa*(eta*P + Q);

    u_c = wrap(mean2true(chief(6),chief(2),tol)+chief(5));

    %% in-plane pair ------------------------------------------------------
    if norm(de)>1e-6 || abs(dl)>1e-6
        xi = atan2(de(1),de(2));
        k  = 0;
        while true
            u1 = xi + k*pi;  u2 = u1 + pi;
            t1 = t_cur + mod(u1-u_c,2*pi)/nbar;
            t2 = t_cur + mod(u2-u_c,2*pi)/nbar;
            if t1>t_cur && t2>t_cur, break; end, k=k+1
        end
        Delta_t   = t2 - t1;
        gamma     = 0.5*J2*(Re/a)^2/eta^4;
        du_J2     = -12*gamma*sin(2*chief(3))*roe_cur(5)*nbar*Delta_t;
        dl_comp   = dl - du_J2;                     % compensated λ

        de_n  = norm(de);
        dvr1 = n*a/2 * (-dl_comp/2 + de_n);
        dvr2 = n*a/2 * (-dl_comp/2 - de_n);
        max_r = 0.5;  dvr1 = sign(dvr1)*min(abs(dvr1),max_r);
                      dvr2 = sign(dvr2)*min(abs(dvr2),max_r);

        burns = [burns;
                 struct('t',t1,'dv_rtn',[dvr1;0;0],'tag','R1','type','reconfig');
                 struct('t',t2,'dv_rtn',[dvr2;0;0],'tag','R2','type','reconfig')];
    end

    %% out-of-plane single pulse -----------------------------------------
    if norm(di) > 1e-6
        if isempty(burns), t_ref = t_cur; else, t_ref = max([burns.t]'); end
        k=0; while true
            uN = atan2(di(2),di(1)) + k*pi;
            tn = t_ref + mod(uN-u_c,2*pi)/nbar;
            if tn>t_ref, break; end, k=k+1
        end
        dvn = n*a*norm(di);   dvn = sign(dvn)*min(abs(dvn),0.5);

        burns = [burns;
                 struct('t',tn,'dv_rtn',[0;0;dvn],'tag','N','type','reconfig')];
    end
end

% -------------------------------------------------------------------------
%  kick, dv_RTN2ECI, get_rtn_components are unchanged from your original
% -------------------------------------------------------------------------
function t_burn = find_next_burn_time(u_tgt, u_now, t_now, nbar)
    % Find the next time when the satellite reaches the target argument of latitude
    
    % Helper functions
    wrap = @(x) wrapTo2Pi(mod(x,2*pi));
    
    k = 0;
    t_burn = -100;
    
    while t_burn < t_now
        u_target = u_tgt + k*2*pi;
        t_burn = t_now + (u_target - u_now)/nbar;
        k = k + 1;
        
        if k > 100
            error('Could not find valid burn time within 100 revolutions');
        end
    end
end

function [chief_new, dep_new] = execute_burn(chief, dep, burn, t_from, dt_fine)
    % Execute a single burn and propagate
    global mu
    
    % Propagate to burn time
    [~, c_tmp] = ode4(@compute_rates_GVE_J2, [t_from burn.t]', chief', dt_fine);
    [~, d_tmp] = ode4(@compute_rates_GVE_J2, [t_from burn.t]', dep', dt_fine);
    chief_at_burn = c_tmp(end,:)';
    dep_at_burn = d_tmp(end,:)';
    
    % Apply burn to deputy
    rv_dep = oe2rv(dep_at_burn', mu);
    r_dep = rv_dep(1:3);
    v_dep = rv_dep(4:6);
    dv_eci = dv_RTN2ECI(r_dep, v_dep, burn.dv_rtn);
    dep_after_burn = kick(dep_at_burn', dv_eci)';
    
    chief_new = chief_at_burn;
    dep_new = dep_after_burn;
end
function oeNew = kick(oeRow, dvECI)
    global mu
    
    % Helper functions
    vrow = @(x) x(:).';        % force row
    
    rv = oe2rv(oeRow, mu);
    rv(4:6) = rv(4:6) + dvECI(:);
    oe = safe_rv2oe(rv, mu);
    oeNew = [oe(1:5) safe_true2mean(oe(6), oe(2))];
end

function rtn_components = get_rtn_components(oe_dep, dv_vec)
    global mu
    
    % Helper functions
    vrow = @(x) x(:).';        % force row
    
    rv_dep = oe2rv(oe_dep, mu);
    r_dep = rv_dep(1:3);
    v_dep = rv_dep(4:6);
    [Rhat, That, Nhat] = eci2rtn_dir(r_dep, v_dep);
    
    dv_R = dot(dv_vec, Rhat);
    dv_T = dot(dv_vec, That);
    dv_N = dot(dv_vec, Nhat);
    
    rtn_components = [dv_R dv_T dv_N];
end

function [R,T,N] = eci2rtn_dir(r, v)
    % Helper functions
    vrow = @(x) x(:).';        % force row
    
    R = r(:).' / norm(r);
    h = cross(r, v);
    N = h(:).' / norm(h);
    T = cross(N, R);
end

function delta_v_ECI = dv_RTN2ECI(r_ECI, v_ECI, delta_v_RTN)
    % Helper functions
    vcol = @(x) x(:);          % force column
    
    % Normalize R = radial direction (along position)
    R_hat = r_ECI / norm(r_ECI);
    R_hat = R_hat(:);

    % N = orbit normal direction = r × v
    N = cross(r_ECI, v_ECI);
    N_hat = N / norm(N);
    N_hat = N_hat(:);

    T_hat = cross(N_hat, R_hat);
    T_hat = T_hat(:);

    % Construct rotation matrix from RTN to ECI
    Q_RTN2ECI = [R_hat, T_hat, N_hat];
    delta_v_ECI = Q_RTN2ECI * delta_v_RTN;
end