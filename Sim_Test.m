function simulate_TDX

clc; clear; close all; format long g

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

%--------------------------------------------------
% 1)  INITIAL CONDITIONS  (chief TerraSAR‑X orbit)
%--------------------------------------------------
a_c    = 6892944.65;                       % m
e_c    = 1.37e-4;
i_c    = deg2rad(97.440124);
RAAN_c = deg2rad(104.274891);
w_c    = deg2rad(67.975723);
M_c    = deg2rad(292.169756);
chief0 = [a_c; e_c; i_c; RAAN_c; w_c; M_c];

%  relative ROE scenario (DEM formation)
rel_qns_init = [0 0 0 300 0 -1000]';     % [m]

dep0   = qns2oe(chief0', rel_qns_init)';   % deputy mean OE

%--------------------------------------------------
% 2)  SIMULATION OPTIONS
%--------------------------------------------------
sim.t0       = 0;                 % s
sim.tf       = 30*86400;          % s  (30 days)
sim.dt       = 60;                % s  integ. step
sim.check_dt = 2*3600;            % s  monitor every 2 h
sim.ei_thr   = 7;                 % deg trigger
sim.goal     = rel_qns_init;      % desired ROE (m)

%--------------------------------------------------
% 3)  PREP: LOGS
%--------------------------------------------------
log.t     = sim.t0;
log.roe   = a_c*compute_roes(chief0',dep0')';
log.burn  = struct('t',{},'dv',{},'roe_post',{});

%--------------------------------------------------
% 4)  MAIN CLOSED‑LOOP
%--------------------------------------------------
t_now      = sim.t0;
t_next_chk = sim.t0 + sim.check_dt;
chief      = chief0;
dep        = dep0;

while t_now < sim.tf

    % ---- propagate until next check (or end) -----------------------
    t_end_seg   = min(t_next_chk, sim.tf);
    [~,chiefSeg]= ode4(@compute_rates_GVE_J2,[t_now t_end_seg]',chief,sim.dt);
    [~,depSeg  ]= ode4(@compute_rates_GVE_J2,[t_now t_end_seg]',dep'  ,sim.dt);
    chief       = chiefSeg(end,:)';
    dep         = depSeg  (end,:)';
    t_now       = t_end_seg;

    % ---- store intermediate samples (for plotting) -----------------
    for k=2:size(chiefSeg,1)
        log.t(end+1,1)= t_now - (size(chiefSeg,1)-k)*sim.dt;
        log.roe(:,end+1)=a_c*compute_roes(chiefSeg(k,:),depSeg(k,:))';
    end

    % ---- at check epochs, test e/i angle ---------------------------
    if abs(t_now - t_next_chk) < 1e-6 && t_now < sim.tf

        roe_now = a_c*compute_roes(chief',dep')';
        ang = ei_angle_deg(roe_now(3:4), roe_now(5:6));

        if ang > sim.ei_thr
            % ~~~ DESIGN TWO‑PULSE MANOEUVRE ~~~
            [dv1T,dv1R,u1, dv2T,dv2R,u2] = ...
                     design_tandem_pair(chief, dep, sim.goal);

            % ---- coast to burn‑1 ----
            [chief,dep,t_b1] = coast_to_u(chief,dep,u1,sim.dt);
            dep = apply_inplane_impulse(dep,dv1T,dv1R);

            % ---- coast to burn‑2 ----
            [chief,dep,~   ] = coast_to_u(chief,dep,u2,sim.dt);
            dep = apply_inplane_impulse(dep,dv2T,dv2R);

            % ---- log burn ----
            log.burn(end+1).t        = t_b1;
            log.burn(end).dv         = [dv1T dv1R; dv2T dv2R];
            log.burn(end).roe_post   = a_c*compute_roes(chief',dep')';
            fprintf('Burn %2d  @ %.1f h  |Δv|=%.2f + %.2f cm/s  angle=%.1f°\n',...
                    numel(log.burn), t_b1/3600, 100*norm([dv1T dv1R]), ...
                    100*norm([dv2T dv2R]), ang);
        end
        t_next_chk = t_next_chk + sim.check_dt;
    end
end

fprintf('\nFinished: %d burns, Σ|Δv| = %.2f cm/s\n', ...
        numel(log.burn), 100*sum(cellfun(@(m)norm(m, 'fro'), {log.burn.dv})));

%------------------------------------------------------------------
% OPTIONAL quick plot
figure; plot(log.t/86400, log.roe(3,:), log.t/86400, log.roe(4,:)); grid on
xlabel('days'); ylabel('a·δe [m]'); legend('\delta e_x','\delta e_y');
title('Deputy e‑vector in closed loop');
end
%% =================================================================
function ang = ei_angle_deg(de, di)
% smallest angle between δe and δi  (deg)
ang = rad2deg(acos( max(-1,min(1,dot(de,di)/(norm(de)*norm(di)))) ));
ang = min(ang, 180-ang);
end
%% =================================================================
function [dvT1,dvR1,u1, dvT2,dvR2,u2] = ...
          design_tandem_pair(chief, dep, roe_goal)
% full TanDEM‑X in‑plane pair (first‑order, antipodal)
global mu
a = chief(1); n = sqrt(mu/a^3);
roe_now = a*compute_roes(chief',dep')';          % current (m)
De  = roe_now(3:4) - roe_goal(3:4);
Da  = roe_now(1)   - roe_goal(1);
u1  = atan2(roe_now(4), roe_now(3));             % burn‑1 at ξ
dvT1= (n*a/4)*( Da + De(1)*cos(u1) - De(2)*sin(u1) );
dvR1=-(n*a/2)*( -Da/2 + De(1)*sin(u1) + De(2)*cos(u1) );
dvT2=-dvT1;  dvR2=-dvR1;                         % antipode
u2 = mod(u1+pi, 2*pi);
end
%% =================================================================
function [chief_out, dep_out, t_out] = ...
          coast_to_u(chief, dep, u_target, dt)
% propagate until chief’s argument of latitude reaches u_target
global mu tol
u_now  = arglat(chief);
du     = mod(u_target - u_now, 2*pi);
n      = sqrt(mu/chief(1)^3);
t_coast= du / n;
[t_seg,chief_traj]=ode4(@compute_rates_GVE_J2,[0 t_coast]',chief,dt);
[~,   dep_traj  ] =ode4(@compute_rates_GVE_J2,[0 t_coast]',dep  ,dt);
chief_out = chief_traj(end,:)';
dep_out   = dep_traj(end,:)';
t_out     = t_seg(end);
end
%% =================================================================
function u = arglat(oe)
% argument of latitude from mean elements
global tol
M = wrapTo2Pi(oe(6)); nu = mean2true(M,oe(2),tol);
u = wrapTo2Pi(nu + oe(5));
end
%% =================================================================
function dep_new = apply_inplane_impulse(dep, dvT, dvR)
% add Δv_T*T̂ + Δv_R*R̂  and convert back to *mean* OE
global mu tol
rv      = oe2rv(dep',mu);
r = rv(1:3); v = rv(4:6);
Rhat = r/norm(r);       What= cross(r,v); That= cross(What,Rhat); That=That/norm(That);
v_new  = v + dvT*That + dvR*Rhat;
oe_true= rv2oe([r;v_new],mu);               % a,e,i,Ω,ω,ν
M_new  = true2mean(oe_true(6), oe_true(2));
dep_new= [oe_true(1:5) M_new]';
end