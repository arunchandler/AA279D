% main.m
clc; clear; close all;
addpath('mean_osc');
format long g;

%% constants
global mu tol
mu  = 3.986004418e14;    % [m^3/s^2]
tol = 1e-10;             % kepler solver tolerance

%% Case A: ecc = 0.1, delta_a = 0
fprintf('\n===== Case A: ecc = 0.1, delta_a = 0 =====\n');
a_caseA    = 6886536.686;     ecc_caseA   = 0.1;
inc_caseA  = deg2rad(97.4453); raan_caseA  = deg2rad(351.0108);
argp_caseA = deg2rad(101.2452); M_caseA    = deg2rad(11.6520);
a_depA     = a_caseA;         ecc_depA    = 0.1 + 5e-6;
inc_depA   = deg2rad(97.4454); raan_depA   = deg2rad(351.0106);
argp_depA  = deg2rad(101.2452);
M_depA     = deg2rad(11.6520);

[rYA_A, rLIN_A, rTR_A, vYA_A, vLIN_A, vTR_A, t_orbit_A] = run_case( ...
  a_caseA,ecc_caseA,inc_caseA,raan_caseA,argp_caseA,M_caseA, ...
  a_depA,ecc_depA,inc_depA,raan_depA,argp_depA,M_depA );

plot_all(...
  rYA_A, rLIN_A, rTR_A, ...
  vYA_A, vLIN_A, vTR_A, ...
  t_orbit_A, 'Case A: ecc=0.1, \delta a=0');

%% Case B: delta_a = –100 m, ecc≈0
fprintf('\n===== Case B: delta_a = -100 m, ecc≈0 =====\n');
a_caseB    = 6886536.686;     ecc_caseB   = 0.1;
inc_caseB  = deg2rad(97.4453); raan_caseB  = deg2rad(351.0108);
argp_caseB = deg2rad(101.2452); M_caseB    = deg2rad(11.6520);
a_depB     = a_caseB - 100;    ecc_depB    = 0.1;
inc_depB   = deg2rad(97.4454); raan_depB   = deg2rad(351.0106);
argp_depB  = deg2rad(-100.5043);
M_depB     = deg2rad(201.086 + 12.35936);

[rYA_B, rLIN_B, rTR_B, vYA_B, vLIN_B, vTR_B, t_orbit_B] = run_case( ...
  a_caseB,ecc_caseB,inc_caseB,raan_caseB,argp_caseB,M_caseB, ...
  a_depB,ecc_depB,inc_depB,raan_depB,argp_depB,M_depB );

plot_all(...
  rYA_B, rLIN_B, rTR_B, ...
  vYA_B, vLIN_B, vTR_B, ...
  t_orbit_B, 'Case B: \delta a=-100 m, ecc≈0');

%% Case C: ecc = 0.6, delta_a = 0
fprintf('\n===== Case C: ecc = 0.6, delta_a = 0 =====\n');
a_caseC    = 6886536.686;     ecc_caseC   = 0.6;
inc_caseC  = deg2rad(97.4453); raan_caseC  = deg2rad(351.0108);
argp_caseC = deg2rad(101.2452); M_caseC    = deg2rad(11.6520);
a_depC     = a_caseC;         ecc_depC    = 0.6;
inc_depC   = deg2rad(97.4454); raan_depC   = deg2rad(351.0106);
argp_depC  = deg2rad(100.5043);
M_depC     = deg2rad(12.35936);

[rYA_C, rLIN_C, rTR_C, vYA_C, vLIN_C, vTR_C, t_orbit_C] = run_case( ...
  a_caseC,ecc_caseC,inc_caseC,raan_caseC,argp_caseC,M_caseC, ...
  a_depC,ecc_depC,inc_depC,raan_depC,argp_depC,M_depC );

plot_all(...
  rYA_C, rLIN_C, rTR_C, ...
  vYA_C, vLIN_C, vTR_C, ...
  t_orbit_C, 'Case C: ecc=0.6, \delta a=0');


%%------------------------------------------------------------------------------%%
function [rYA, rLIN, rTR, vYA, vLIN, vTR, t_orbit] = run_case(...
  a_chief, ecc_chief, inc_chief, raan_chief, argp_chief, M_chief, ...
  a_deputy, ecc_deputy, inc_deputy, raan_deputy, argp_deputy, M_deputy)

  global mu tol

  %% a) initial states & print
  nu0    = mean2true(M_chief, ecc_chief, tol);
  rvC    = oe2rv([a_chief,ecc_chief,inc_chief,raan_chief,argp_chief,M_chief], mu);
  rvD    = oe2rv([a_deputy,ecc_deputy,inc_deputy,raan_deputy,argp_deputy,M_deputy], mu);
  [rv_rel0,~] = eci2rtn(rvC, rvD);

  a_diff = a_deputy - a_chief;
  inc_diff = inc_deputy - inc_chief;
  ecc_diff = ecc_deputy - ecc_chief;
  raan_diff = raan_deputy - raan_chief;
  argp_diff = argp_deputy - argp_chief;
  M_diff = M_deputy - M_chief;

  r_peri = a_chief*(1 - ecc_chief);
  ro_init = norm(rv_rel0(1:3));
  p       = a_chief*(1-ecc_chief^2);
  r0      = p/(1+ecc_chief*cos(nu0));
  ratio0  = ro_init/r_peri;

  fprintf('\nRelevant Initial Orbital Elements for TH Formulation\n');
  fprintf('---------------------------------------------------\n');
  fprintf('d a           : %+12.6f m\n', a_diff);
  fprintf('d i           : %+12.6f deg\n', rad2deg(inc_diff));
  fprintf('d e           : %+12.6e\n', ecc_diff);
  fprintf('d RAAN        : %+12.6f deg\n', rad2deg(raan_diff));
  fprintf('d argp        : %+12.6f deg\n', rad2deg(argp_diff));
  fprintf('d M           : %+12.6f deg\n', rad2deg(M_diff));
  fprintf('Initial r_peri    : %12.6f m\n', r_peri);
  fprintf('Initial ro_init   : %12.6f m\n', ro_init);
  fprintf('Initial r0        : %12.6f m\n', r0);
  fprintf('Initial ratio ro/r_peri : %.3e\n', ratio0);
  fprintf('Chief ecc         : %.6f\n', ecc_chief);

  %% b) YA constants & print
  Ks = getYAconstants(rv_rel0, a_chief, ecc_chief, nu0);
  %Ks(1) = 0;
  fprintf('\nYA Integration Constants:\n');
  fprintf('-------------------------\n');
  for ii = 1:6
      fprintf('K%d = %12.6e\n', ii, Ks(ii));
  end
  
  %% c) time grid
  n   = sqrt(mu/a_chief^3);
  T   = 2*pi/n;
  N   = 10000;
  tvec = linspace(0,15*T,N).';
  t_orbit = tvec/T;
  eta = sqrt(1-ecc_chief^2);

  %% d) YA analytic
  rYA = zeros(N,3); vYA = zeros(N,3);
  for k=1:N
    t    = tvec(k);
    M    = M_chief + n*t;
    f    = mean2true(M, ecc_chief, tol);
    tau  = n*t/eta^3;
    Phi  = buildYAphi(a_chief,ecc_chief,f, tau);
    X    = Phi * Ks;
    rYA(k,:) = X(1:3)';
    vYA(k,:) = X(4:6)';
  end

  %% e) geometric linear map
  da = (a_deputy-a_chief)/a_chief;
  dl = (M_deputy+argp_deputy - M_chief - argp_chief) + (raan_deputy-raan_chief)*cos(inc_chief);
  ex = ecc_deputy*cos(argp_deputy) - ecc_chief*cos(argp_chief);
  ey = ecc_deputy*sin(argp_deputy) - ecc_chief*sin(argp_chief);
  ix = inc_deputy-inc_chief;
  iy = (raan_deputy-raan_chief)*sin(inc_chief);
  ROE = [da;dl;ex;ey;ix;iy];

  %--- print ROE ---
  fprintf('\nRelative Quasi-Nonsingular Orbit Elements:\n');
  fprintf('---------------------------------------------------\n');
  fprintf('da/a   : %12.6e\n', da);
  fprintf('dλ     : %12.6e\n', dl);
  fprintf('de_x   : %12.6e\n', ex);
  fprintf('de_y   : %12.6e\n', ey);
  fprintf('di_x   : %12.6e\n', ix);
  fprintf('di_y   : %12.6e\n', iy);
  
  [rv_LIN] = propagateLinearEcc(ROE, ...
      [a_chief; ecc_chief; inc_chief; raan_chief; argp_chief; M_chief], ...
      tvec);

  rLIN = rv_LIN(:, 1:3);
  vLIN = rv_LIN(:, 4:6); 


  %% f) “truth” ODE & print
  dt    = 15*T/(N-1); 
  [tout, S] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [0;15*T], [rv_rel0;rvC], dt);
  rTR = S(:,1:3); vTR = S(:,4:6);

  rho = vecnorm(rTR,2,2);
  rChief = zeros(N,1);
  for k=1:N
    Mk = M_chief + n*tout(k);
    Ek = mean2ecc(Mk,ecc_chief,tol);
    fk = ecc2true(Ek,ecc_chief);
    rChief(k) = a_chief*(1-ecc_chief^2)/(1+ecc_chief*cos(fk));
  end
  ratio = rho./rChief;
  fprintf('\nNonlinear Trajectory Separation Check\n');
  fprintf('-------------------------------------\n');
  fprintf('Mean ||rho||/rChief = %.3e\n', mean(ratio));
  fprintf('Max  ||rho||/rChief = %.3e\n', max(ratio));
  fprintf('Min  ||rho||/rChief = %.3e\n', min(ratio));

  % convert to km / km-s
  %rYA  = rYA/1e3;  vYA  = vYA/1e3;
  %rLIN = rLIN/1e3; vLIN = vLIN/1e3;
  %rTR  = rTR/1e3;  vTR  = vTR/1e3;
end

function plot_all(rYA,rLIN,rTR, vYA,vLIN,vTR, t_orbit, titleStr)

  %% 1) all three overlay—planes
  figure('Name',[titleStr '—All 3'],'NumberTitle','off');
  subplot(1,3,1)
    plot(rYA(:,2),rYA(:,1),'r-','LineWidth',1.0); hold on;
    plot(rLIN(:,2),rLIN(:,1),'b:','LineWidth',1.0);
    plot(rTR(:,2),rTR(:,1),'c:','LineWidth',1.0);
    grid on;  xlabel('T [m]'); ylabel('R [m]');
    legend('YA','Geo Map','Diff eq','Location','best');
  subplot(1,3,2)
    plot(rYA(:,3),rYA(:,1),'r-','LineWidth',1.0); hold on;
    plot(rLIN(:,3),rLIN(:,1),'b:','LineWidth',1.0);
    plot(rTR(:,3),rTR(:,1),'c:','LineWidth',1.0);
    grid on;  xlabel('N [m]'); ylabel('R [m]');
  subplot(1,3,3)
    plot(rYA(:,2),rYA(:,3),'r-','LineWidth',1.0); hold on;
    plot(rLIN(:,2),rLIN(:,3),'b:','LineWidth',1.0);
    plot(rTR(:,2),rTR(:,3),'c:','LineWidth',1.0);
    grid on;  xlabel('T [m]'); ylabel('N [m]');

  %% 2) all three overlay—3D
  figure('Name',[titleStr '—All 3 (3D)'],'NumberTitle','off');
  plot3(rYA(:,1),rYA(:,2),rYA(:,3),'r-','LineWidth',1.0); hold on;
  plot3(rLIN(:,1),rLIN(:,2),rLIN(:,3),'b:','LineWidth',1.0);
  plot3(rTR(:,1),rTR(:,2),rTR(:,3),'c:','LineWidth',1.0);
  grid on; view(3);
  xlabel('R [m]'); ylabel('T [m]'); zlabel('N [m]');
  legend('YA','Geo Map','Diff eq','Location','best');

  %% 3) velocity vs orbits
  figure('Name',[titleStr '—Vel vs Orbits'],'NumberTitle','off');
  subplot(3,1,1)
    plot(t_orbit,vYA(:,1),'r-','LineWidth',1.0); hold on;
    plot(t_orbit,vLIN(:,1),'b:','LineWidth',1.0);
    plot(t_orbit,vTR(:,1),'c:','LineWidth',1.0);
    grid on; xlabel('Orbits'); ylabel('Ṙ [m/s]');
    legend('YA','Geo Map','Diff eq','Location','best');
  subplot(3,1,2)
    plot(t_orbit,vYA(:,2),'r-','LineWidth',1.0); hold on;
    plot(t_orbit,vLIN(:,2),'b:','LineWidth',1.0);
    plot(t_orbit,vTR(:,2),'c:','LineWidth',1.0);
    grid on; xlabel('Orbits'); ylabel('Ṫ [m/s]');
  subplot(3,1,3)
    plot(t_orbit,vYA(:,3),'r-','LineWidth',1.0); hold on;
    plot(t_orbit,vLIN(:,3),'b:','LineWidth',1.0);
    plot(t_orbit,vTR(:,3),'c:','LineWidth',1.0);
    grid on; xlabel('Orbits'); ylabel('Ṅ [m/s]');

  %% 4) propagation errors vs orbits
  errYA  = rYA - rTR;
  errLIN = rLIN - rTR;
  figure('Name',[titleStr '—Prop Errors'],'NumberTitle','off');
  subplot(3,1,1)
    plot(t_orbit,errYA(:,1),'r-','LineWidth',1.0); hold on;
    plot(t_orbit,errLIN(:,1),'b:','LineWidth',1.0);
    grid on; xlabel('Orbits'); ylabel('R error [m]');
    legend('YA–Truth','Geo–Truth','Location','best');
  subplot(3,1,2)
    plot(t_orbit,errYA(:,2),'r-','LineWidth',1.0); hold on;
    plot(t_orbit,errLIN(:,2),'b:','LineWidth',1.0);
    grid on; xlabel('Orbits'); ylabel('T error [m]');
  subplot(3,1,3)
    plot(t_orbit,errYA(:,3),'r-','LineWidth',1.0); hold on;
    plot(t_orbit,errLIN(:,3),'b:','LineWidth',1.0);
    grid on; xlabel('Orbits'); ylabel('N error [m]');
end
