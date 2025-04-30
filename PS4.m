clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

% Part a) initial conditions

%initial cheif elements
a_TSX_init_1 = 6886536.686;  % m
e_TSX_init_1 = 0.0001264;
i_TSX_init_1 = deg2rad(97.4453);
RAAN_TSX_init_1 = deg2rad(351.0108);
omega_TSX_init_1 = deg2rad(101.2452);
M_TSX_init_1 = deg2rad(11.6520);
nu_TSX_init_1 = mean2true(M_TSX_init_1, e_TSX_init_1,tol);
u_TSX_init_1 = nu_TSX_init_1 + omega_TSX_init_1;

%initial deputy elements - changed to have similar separation for R,T,N
a_TDX_init_1 =  6886536.686;
e_TDX_init_1 = 0.0001269;
i_TDX_init_1 = deg2rad(97.4554);
RAAN_TDX_init_1 = deg2rad(351.0106);
omega_TDX_init_1 = deg2rad(-100.5043);
M_TDX_init_1 = deg2rad(201.086+12.34926);

%timing parameters
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init_1^3);
n_orbit    = 15;
T          = 2*pi / n;
tend       = n_orbit*T;
num_points = 1000;
dt         = (tend - tstart)/(num_points-1);
t_grid     = linspace(tstart, tend, num_points).';    % [s]
t_orbit    = t_grid / T;

TSX_init_oe_1 = [a_TSX_init_1, e_TSX_init_1, i_TSX_init_1, RAAN_TSX_init_1, omega_TSX_init_1, M_TSX_init_1];
TSX_init_osc_oe_1 = TSX_init_oe_1;

TDX_init_oe_1 = [a_TDX_init_1, e_TDX_init_1, i_TDX_init_1, RAAN_TDX_init_1, omega_TDX_init_1, M_TDX_init_1];
TDX_init_osc_oe_1 = TDX_init_oe_1;

TSX_init_rv_1 = oe2rv(TSX_init_osc_oe_1, mu)
TDX_init_rv_1 = oe2rv(TDX_init_osc_oe_1, mu)

[rel_rv_1,~] = eci2rtn(TSX_init_rv_1, TDX_init_rv_1)

%Part b) Initial condiditions from quasi-nonsingular

%initial cheif elements & state
a_TSX_init_2 = a_TSX_init_1;
e_TSX_init_2 = e_TSX_init_1;
i_TSX_init_2 = i_TSX_init_1;
RAAN_TSX_init_2 = RAAN_TSX_init_1;
omega_TSX_init_2 = omega_TSX_init_1;
M_TSX_init_2 = M_TSX_init_1;
nu_TSX_init_2 = nu_TSX_init_1;
u_TSX_init_2 = u_TSX_init_1;

TSX_init_oe_2 = [a_TSX_init_2, e_TSX_init_2, i_TSX_init_2, RAAN_TSX_init_2, omega_TSX_init_2, M_TSX_init_2];
TSX_init_osc_oe_2 = TSX_init_oe_2;

TSX_init_rv_2 = oe2rv(TSX_init_osc_oe_2, mu)

%quasi-nonsingular elements
rel_qns = [0, 100, 50, 100, 30, 200];

%compute deputy elements from quasi-nonsingular
TDX_init_oe_2 = qns2oe(TSX_init_oe_2, rel_qns)
TDX_init_osc_oe_2 = TDX_init_oe_2;
TDX_init_rv_2 = oe2rv(TDX_init_osc_oe_2, mu)

[rel_rv_2,~] = eci2rtn(TSX_init_rv_2, TDX_init_rv_2)

% Part c) numerical integration

[t_out, TSX_rv_1]    = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', TSX_init_rv_1, dt);
[t_out, TDX_rv_1]    = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', TDX_init_rv_1, dt);
[t_out, TSX_rv_J2_1] = ode4(@compute_rates_rv_perturbed, [tstart, tend]', TSX_init_rv_1, dt);
[t_out, TDX_rv_J2_1] = ode4(@compute_rates_rv_perturbed, [tstart, tend]', TDX_init_rv_1, dt);
[t_out, TSX_rv_2]    = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', TSX_init_rv_2, dt);
[t_out, TDX_rv_2]    = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', TDX_init_rv_2, dt);
[t_out, TSX_rv_J2_2] = ode4(@compute_rates_rv_perturbed, [tstart, tend]', TSX_init_rv_2, dt);
[t_out, TDX_rv_J2_2] = ode4(@compute_rates_rv_perturbed, [tstart, tend]', TDX_init_rv_2, dt);

TSX_osc_1       = zeros(num_points, 6); %[a, ex, ey, i, RAAN, u]
TSX_mean_1      = zeros(num_points, 6); %[a, ex, ey, i, RAAN, u]
TDX_osc_1       = zeros(num_points, 6);
TDX_mean_1      = zeros(num_points, 6);
rel_osc_qns_1   = zeros(num_points, 6); %[da, dlambda, dex, dey, dix, diy]
rel_mean_qns_1  = zeros(num_points, 6); %[da, dlambda, dex, dey, dix, diy]

TSX_J2_osc_1    = zeros(num_points, 6);
TSX_J2_mean_1   = zeros(num_points, 6);
TDX_J2_osc_1    = zeros(num_points, 6);
TDX_J2_mean_1   = zeros(num_points, 6);
rel_J2_osc_qns_1= zeros(num_points, 6);
rel_J2_mean_qns_1= zeros(num_points, 6);

TSX_osc_2       = zeros(num_points, 6);
TSX_mean_2      = zeros(num_points, 6);
TDX_osc_2       = zeros(num_points, 6);
TDX_mean_2      = zeros(num_points, 6);
rel_osc_qns_2   = zeros(num_points, 6);
rel_mean_qns_2  = zeros(num_points, 6);

TSX_J2_osc_2    = zeros(num_points, 6);
TSX_J2_mean_2   = zeros(num_points, 6);
TDX_J2_osc_2    = zeros(num_points, 6);
TDX_J2_mean_2   = zeros(num_points, 6);
rel_J2_osc_qns_2= zeros(num_points, 6);
rel_J2_mean_qns_2= zeros(num_points, 6);

TSX_params_1      = zeros(12);
TDX_params_1      = zeros(12);
TSX_params_J2_1   = zeros(12);
TDX_params_J2_1   = zeros(12);
TSX_params_2      = zeros(12);
TDX_params_2      = zeros(12);
TSX_params_J2_2   = zeros(12);
TDX_params_J2_2   = zeros(12);

for idx = 1:num_points

    % --- Set 1, no J2 ---
    TSX_params_1      = rv2oe(TSX_rv_1(idx,:),    mu);
    TDX_params_1      = rv2oe(TDX_rv_1(idx,:),    mu);

    a_TSX_osc_1       = TSX_params_1(1);
    e_TSX_osc_1       = TSX_params_1(2);
    i_TSX_osc_1       = TSX_params_1(3);
    RAAN_TSX_osc_1    = TSX_params_1(4);
    omega_TSX_osc_1   = TSX_params_1(5);
    nu_TSX_osc_1      = TSX_params_1(6);
    e_vec_TSX_osc_1   = TSX_params_1(11:12);
    M_TSX_osc_1       = true2mean(nu_TSX_osc_1, e_TSX_osc_1);
    u_TSX_osc_1       = M_TSX_osc_1 + omega_TSX_osc_1;
    TSX_osc_1(idx,:)  = [a_TSX_osc_1, e_vec_TSX_osc_1(:)', i_TSX_osc_1, RAAN_TSX_osc_1, u_TSX_osc_1];

    TSX_mean_oe_1     = osc2mean([a_TSX_osc_1, e_TSX_osc_1, i_TSX_osc_1, ...
                                 RAAN_TSX_osc_1, omega_TSX_osc_1, M_TSX_osc_1],0)'; %0 for no J2
    a_TSX_mean_1      = TSX_mean_oe_1(1);
    e_TSX_mean_1      = TSX_mean_oe_1(2);
    i_TSX_mean_1      = TSX_mean_oe_1(3);
    RAAN_TSX_mean_1   = wrapTo2Pi(TSX_mean_oe_1(4));
    omega_TSX_mean_1  = wrapTo2Pi(TSX_mean_oe_1(5));
    M_TSX_mean_1      = TSX_mean_oe_1(6);
    ex_TSX_mean_1     = e_TSX_mean_1*cos(omega_TSX_mean_1);
    ey_TSX_mean_1     = e_TSX_mean_1*sin(omega_TSX_mean_1);
    u_TSX_mean_1      = M_TSX_mean_1 + omega_TSX_mean_1;
    TSX_mean_1(idx,:) = [a_TSX_mean_1, ex_TSX_mean_1, ey_TSX_mean_1, ...
                         i_TSX_mean_1, RAAN_TSX_mean_1, u_TSX_mean_1];

    a_TDX_osc_1       = TDX_params_1(1);
    e_TDX_osc_1       = TDX_params_1(2);
    i_TDX_osc_1       = TDX_params_1(3);
    RAAN_TDX_osc_1    = TDX_params_1(4);
    omega_TDX_osc_1   = TDX_params_1(5);
    nu_TDX_osc_1      = TDX_params_1(6);
    e_vec_TDX_osc_1   = TDX_params_1(11:12);
    M_TDX_osc_1       = true2mean(nu_TDX_osc_1, e_TDX_osc_1);
    u_TDX_osc_1       = M_TDX_osc_1 + omega_TDX_osc_1;
    TDX_osc_1(idx,:)  = [a_TDX_osc_1, e_vec_TDX_osc_1(:)', i_TDX_osc_1, ...
                         RAAN_TDX_osc_1, u_TDX_osc_1];

    TDX_mean_oe_1     = osc2mean([a_TDX_osc_1, e_TDX_osc_1, i_TDX_osc_1, ...
                                 RAAN_TDX_osc_1, omega_TDX_osc_1, M_TDX_osc_1],0)';%0 for no J2
    a_TDX_mean_1      = TDX_mean_oe_1(1);
    e_TDX_mean_1      = TDX_mean_oe_1(2);
    i_TDX_mean_1      = TDX_mean_oe_1(3);
    RAAN_TDX_mean_1   = wrapTo2Pi(TDX_mean_oe_1(4));
    omega_TDX_mean_1  = wrapTo2Pi(TDX_mean_oe_1(5));
    M_TDX_mean_1      = TDX_mean_oe_1(6);
    ex_TDX_mean_1     = e_TDX_mean_1*cos(omega_TDX_mean_1);
    ey_TDX_mean_1     = e_TDX_mean_1*sin(omega_TDX_mean_1);
    u_TDX_mean_1      = M_TDX_mean_1 + omega_TDX_mean_1;
    TDX_mean_1(idx,:) = [a_TDX_mean_1, ex_TDX_mean_1, ey_TDX_mean_1, ...
                         i_TDX_mean_1, RAAN_TDX_mean_1, u_TDX_mean_1];

    rel_osc_qns_1(idx,:) = a_TSX_osc_1*compute_roes(a_TSX_osc_1, i_TSX_osc_1, e_TSX_osc_1, ...
                                        RAAN_TSX_osc_1, u_TSX_osc_1, ...
                                        a_TDX_osc_1, i_TDX_osc_1, e_TDX_osc_1, ...
                                        RAAN_TDX_osc_1, u_TDX_osc_1);

    rel_mean_qns_1(idx,:)= a_TSX_mean_1*compute_roes(a_TSX_mean_1, i_TSX_mean_1, e_TSX_mean_1, ...
                                        RAAN_TSX_mean_1, u_TSX_mean_1, ...
                                        a_TDX_mean_1, i_TDX_mean_1, e_TDX_mean_1, ...
                                        RAAN_TDX_mean_1, u_TDX_mean_1);


    % --- Set 1, with J2 ---
    TSX_params_J2_1   = rv2oe(TSX_rv_J2_1(idx,:), mu);
    TDX_params_J2_1   = rv2oe(TDX_rv_J2_1(idx,:), mu);

    a_TSX_J2_1       = TSX_params_J2_1(1);
    e_TSX_J2_1       = TSX_params_J2_1(2);
    i_TSX_J2_1       = TSX_params_J2_1(3);
    RAAN_TSX_J2_1    = TSX_params_J2_1(4);
    omega_TSX_J2_1   = TSX_params_J2_1(5);
    nu_TSX_J2_1      = TSX_params_J2_1(6);
    e_vec_TSX_J2_1   = TSX_params_J2_1(11:12);
    M_TSX_J2_1       = true2mean(nu_TSX_J2_1, e_TSX_J2_1);
    u_TSX_J2_1       = M_TSX_J2_1 + omega_TSX_J2_1;
    TSX_J2_osc_1(idx,:)= [a_TSX_J2_1, e_vec_TSX_J2_1(:)', i_TSX_J2_1, ...
                          RAAN_TSX_J2_1, u_TSX_J2_1];

    TSX_J2_mean_oe_1  = osc2mean([a_TSX_J2_1, e_TSX_J2_1, i_TSX_J2_1, ...
                                 RAAN_TSX_J2_1, omega_TSX_J2_1, M_TSX_J2_1],1)'; %1 for J2
    a_TSX_J2_mean_1   = TSX_J2_mean_oe_1(1);
    e_TSX_J2_mean_1   = TSX_J2_mean_oe_1(2);
    i_TSX_J2_mean_1   = TSX_J2_mean_oe_1(3);
    RAAN_TSX_J2_mean_1= wrapTo2Pi(TSX_J2_mean_oe_1(4));
    omega_TSX_J2_mean_1=wrapTo2Pi(TSX_J2_mean_oe_1(5));
    M_TSX_J2_mean_1   = TSX_J2_mean_oe_1(6);
    ex_TSX_J2_mean_1  = e_TSX_J2_mean_1*cos(omega_TSX_J2_mean_1);
    ey_TSX_J2_mean_1  = e_TSX_J2_mean_1*sin(omega_TSX_J2_mean_1);
    u_TSX_J2_mean_1   = M_TSX_J2_mean_1 + omega_TSX_J2_mean_1;
    TSX_J2_mean_1(idx,:)= [a_TSX_J2_mean_1, ex_TSX_J2_mean_1, ey_TSX_J2_mean_1, ...
                           i_TSX_J2_mean_1, RAAN_TSX_J2_mean_1, u_TSX_J2_mean_1];

    a_TDX_J2_1       = TDX_params_J2_1(1);
    e_TDX_J2_1       = TDX_params_J2_1(2);
    i_TDX_J2_1       = TDX_params_J2_1(3);
    RAAN_TDX_J2_1    = TDX_params_J2_1(4);
    omega_TDX_J2_1   = TDX_params_J2_1(5);
    nu_TDX_J2_1      = TDX_params_J2_1(6);
    e_vec_TDX_J2_1   = TDX_params_J2_1(11:12);
    M_TDX_J2_1       = true2mean(nu_TDX_J2_1, e_TDX_J2_1);
    u_TDX_J2_1       = M_TDX_J2_1 + omega_TDX_J2_1;
    TDX_J2_osc_1(idx,:)= [a_TDX_J2_1, e_vec_TDX_J2_1(:)', i_TDX_J2_1, ...
                          RAAN_TDX_J2_1, u_TDX_J2_1];

    TDX_J2_mean_oe_1  = osc2mean([a_TDX_J2_1, e_TDX_J2_1, i_TDX_J2_1, ...
                                 RAAN_TDX_J2_1, omega_TDX_J2_1, M_TDX_J2_1],1)';%1 for J2
    a_TDX_J2_mean_1   = TDX_J2_mean_oe_1(1);
    e_TDX_J2_mean_1   = TDX_J2_mean_oe_1(2);
    i_TDX_J2_mean_1   = TDX_J2_mean_oe_1(3);
    RAAN_TDX_J2_mean_1= wrapTo2Pi(TDX_J2_mean_oe_1(4));
    omega_TDX_J2_mean_1=wrapTo2Pi(TDX_J2_mean_oe_1(5));
    M_TDX_J2_mean_1   = TDX_J2_mean_oe_1(6);
    ex_TDX_J2_mean_1  = e_TDX_J2_mean_1*cos(omega_TDX_J2_mean_1);
    ey_TDX_J2_mean_1  = e_TDX_J2_mean_1*sin(omega_TDX_J2_mean_1);
    u_TDX_J2_mean_1   = M_TDX_J2_mean_1 + omega_TDX_J2_mean_1;
    TDX_J2_mean_1(idx,:)= [a_TDX_J2_mean_1, ex_TDX_J2_mean_1, ey_TDX_J2_mean_1, ...
                           i_TDX_J2_mean_1, RAAN_TDX_J2_mean_1, u_TDX_J2_mean_1];

    rel_J2_osc_qns_1(idx,:)= a_TSX_J2_1*compute_roes(a_TSX_J2_1, i_TSX_J2_1, e_TSX_J2_1, ...
                                          RAAN_TSX_J2_1, u_TSX_J2_1, ...
                                          a_TDX_J2_1, i_TDX_J2_1, e_TDX_J2_1, ...
                                          RAAN_TDX_J2_1, u_TDX_J2_1);
    rel_J2_mean_qns_1(idx,:)=a_TSX_J2_mean_1*compute_roes(a_TSX_J2_mean_1, i_TSX_J2_mean_1, e_TSX_J2_mean_1, ...
                                          RAAN_TSX_J2_mean_1, u_TSX_J2_mean_1, ...
                                          a_TDX_J2_mean_1, i_TDX_J2_mean_1, e_TDX_J2_mean_1, ...
                                          RAAN_TDX_J2_mean_1, u_TDX_J2_mean_1);


    % --- Set 2, no J2 ---
    TSX_params_2      = rv2oe(TSX_rv_2(idx,:),    mu);
    TDX_params_2      = rv2oe(TDX_rv_2(idx,:),    mu);

    a_TSX_osc_2       = TSX_params_2(1);
    e_TSX_osc_2       = TSX_params_2(2);
    i_TSX_osc_2       = TSX_params_2(3);
    RAAN_TSX_osc_2    = TSX_params_2(4);
    omega_TSX_osc_2   = TSX_params_2(5);
    nu_TSX_osc_2      = TSX_params_2(6);
    e_vec_TSX_osc_2   = TSX_params_2(11:12);
    M_TSX_osc_2       = true2mean(nu_TSX_osc_2, e_TSX_osc_2);
    u_TSX_osc_2       = M_TSX_osc_2 + omega_TSX_osc_2;
    TSX_osc_2(idx,:)  = [a_TSX_osc_2, e_vec_TSX_osc_2(:)', i_TSX_osc_2, ...
                         RAAN_TSX_osc_2, u_TSX_osc_2];

    TSX_mean_oe_2     = osc2mean([a_TSX_osc_2, e_TSX_osc_2, i_TSX_osc_2, ...
                                 RAAN_TSX_osc_2, omega_TSX_osc_2, M_TSX_osc_2],0)';
    a_TSX_mean_2      = TSX_mean_oe_2(1);
    e_TSX_mean_2      = TSX_mean_oe_2(2);
    i_TSX_mean_2      = TSX_mean_oe_2(3);
    RAAN_TSX_mean_2   = wrapTo2Pi(TSX_mean_oe_2(4));
    omega_TSX_mean_2  = wrapTo2Pi(TSX_mean_oe_2(5));
    M_TSX_mean_2      = TSX_mean_oe_2(6);
    ex_TSX_mean_2     = e_TSX_mean_2*cos(omega_TSX_mean_2);
    ey_TSX_mean_2     = e_TSX_mean_2*sin(omega_TSX_mean_2);
    u_TSX_mean_2      = M_TSX_mean_2 + omega_TSX_mean_2;
    TSX_mean_2(idx,:) = [a_TSX_mean_2, ex_TSX_mean_2, ey_TSX_mean_2, ...
                         i_TSX_mean_2, RAAN_TSX_mean_2, u_TSX_mean_2];

    a_TDX_osc_2       = TDX_params_2(1);
    e_TDX_osc_2       = TDX_params_2(2);
    i_TDX_osc_2       = TDX_params_2(3);
    RAAN_TDX_osc_2    = TDX_params_2(4);
    omega_TDX_osc_2   = TDX_params_2(5);
    nu_TDX_osc_2      = TDX_params_2(6);
    e_vec_TDX_osc_2   = TDX_params_2(11:12);
    M_TDX_osc_2       = true2mean(nu_TDX_osc_2, e_TDX_osc_2);
    u_TDX_osc_2       = M_TDX_osc_2 + omega_TDX_osc_2;
    TDX_osc_2(idx,:)  = [a_TDX_osc_2, e_vec_TDX_osc_2(:)', i_TDX_osc_2, ...
                         RAAN_TDX_osc_2, u_TDX_osc_2];

    TDX_mean_oe_2     = osc2mean([a_TDX_osc_2, e_TDX_osc_2, i_TDX_osc_2, ...
                                 RAAN_TDX_osc_2, omega_TDX_osc_2, M_TDX_osc_2],0)';
    a_TDX_mean_2      = TDX_mean_oe_2(1);
    e_TDX_mean_2      = TDX_mean_oe_2(2);
    i_TDX_mean_2      = TDX_mean_oe_2(3);
    RAAN_TDX_mean_2   = wrapTo2Pi(TDX_mean_oe_2(4));
    omega_TDX_mean_2  = wrapTo2Pi(TDX_mean_oe_2(5));
    M_TDX_mean_2      = TDX_mean_oe_2(6);
    ex_TDX_mean_2     = e_TDX_mean_2*cos(omega_TDX_mean_2);
    ey_TDX_mean_2     = e_TDX_mean_2*sin(omega_TDX_mean_2);
    u_TDX_mean_2      = M_TDX_mean_2 + omega_TDX_mean_2;
    TDX_mean_2(idx,:) = [a_TDX_mean_2, ex_TDX_mean_2, ey_TDX_mean_2, ...
                         i_TDX_mean_2, RAAN_TDX_mean_2, u_TDX_mean_2];

    rel_osc_qns_2(idx,:)= a_TSX_osc_2*compute_roes(a_TSX_osc_2, i_TSX_osc_2, e_TSX_osc_2, ...
                                        RAAN_TSX_osc_2, u_TSX_osc_2, ...
                                        a_TDX_osc_2, i_TDX_osc_2, e_TDX_osc_2, ...
                                        RAAN_TDX_osc_2, u_TDX_osc_2);
    rel_mean_qns_2(idx,:)=a_TSX_mean_2*compute_roes(a_TSX_mean_2, i_TSX_mean_2, e_TSX_mean_2, ...
                                       RAAN_TSX_mean_2, u_TSX_mean_2, ...
                                       a_TDX_mean_2, i_TDX_mean_2, e_TDX_mean_2, ...
                                       RAAN_TDX_mean_2, u_TDX_mean_2);


    % --- Set 2, with J2 ---
    TSX_params_J2_2   = rv2oe(TSX_rv_J2_2(idx,:), mu);
    TDX_params_J2_2   = rv2oe(TDX_rv_J2_2(idx,:), mu);

    a_TSX_J2_2       = TSX_params_J2_2(1);
    e_TSX_J2_2       = TSX_params_J2_2(2);
    i_TSX_J2_2       = TSX_params_J2_2(3);
    RAAN_TSX_J2_2    = TSX_params_J2_2(4);
    omega_TSX_J2_2   = TSX_params_J2_2(5);
    nu_TSX_J2_2      = TSX_params_J2_2(6);
    e_vec_TSX_J2_2   = TSX_params_J2_2(11:12);
    M_TSX_J2_2       = true2mean(nu_TSX_J2_2, e_TSX_J2_2);
    u_TSX_J2_2       = M_TSX_J2_2 + omega_TSX_J2_2;
    TSX_J2_osc_2(idx,:)= [a_TSX_J2_2, e_vec_TSX_J2_2(:)', i_TSX_J2_2, ...
                          RAAN_TSX_J2_2, u_TSX_J2_2];

    TSX_J2_mean_oe_2  = osc2mean([a_TSX_J2_2, e_TSX_J2_2, i_TSX_J2_2, ...
                                 RAAN_TSX_J2_2, omega_TSX_J2_2, M_TSX_J2_2],1)'; %1 for J2
    a_TSX_J2_mean_2   = TSX_J2_mean_oe_2(1);
    e_TSX_J2_mean_2   = TSX_J2_mean_oe_2(2);
    i_TSX_J2_mean_2   = TSX_J2_mean_oe_2(3);
    RAAN_TSX_J2_mean_2= wrapTo2Pi(TSX_J2_mean_oe_2(4));
    omega_TSX_J2_mean_2=wrapTo2Pi(TSX_J2_mean_oe_2(5));
    M_TSX_J2_mean_2   = TSX_J2_mean_oe_2(6);
    ex_TSX_J2_mean_2  = e_TSX_J2_mean_2*cos(omega_TSX_J2_mean_2);
    ey_TSX_J2_mean_2  = e_TSX_J2_mean_2*sin(omega_TSX_J2_mean_2);
    u_TSX_J2_mean_2   = M_TSX_J2_mean_2 + omega_TSX_J2_mean_2;
    TSX_J2_mean_2(idx,:)= [a_TSX_J2_mean_2, ex_TSX_J2_mean_2, ey_TSX_J2_mean_2, ...
                           i_TSX_J2_mean_2, RAAN_TSX_J2_mean_2, u_TSX_J2_mean_2];

    a_TDX_J2_2       = TDX_params_J2_2(1);
    e_TDX_J2_2       = TDX_params_J2_2(2);
    i_TDX_J2_2       = TDX_params_J2_2(3);
    RAAN_TDX_J2_2    = TDX_params_J2_2(4);
    omega_TDX_J2_2   = TDX_params_J2_2(5);
    nu_TDX_J2_2      = TDX_params_J2_2(6);
    e_vec_TDX_J2_2   = TDX_params_J2_2(11:12);
    M_TDX_J2_2       = true2mean(nu_TDX_J2_2, e_TDX_J2_2);
    u_TDX_J2_2       = M_TDX_J2_2 + omega_TDX_J2_2;
    TDX_J2_osc_2(idx,:)= [a_TDX_J2_2, e_vec_TDX_J2_2(:)', i_TDX_J2_2, ...
                          RAAN_TDX_J2_2, u_TDX_J2_2];

    TDX_J2_mean_oe_2  = osc2mean([a_TDX_J2_2, e_TDX_J2_2, i_TDX_J2_2, ...
                                 RAAN_TDX_J2_2, omega_TDX_J2_2, M_TDX_J2_2],1)'; %1 for J2
    a_TDX_J2_mean_2   = TDX_J2_mean_oe_2(1);
    e_TDX_J2_mean_2   = TDX_J2_mean_oe_2(2);
    i_TDX_J2_mean_2   = TDX_J2_mean_oe_2(3);
    RAAN_TDX_J2_mean_2= wrapTo2Pi(TDX_J2_mean_oe_2(4));
    omega_TDX_J2_mean_2=wrapTo2Pi(TDX_J2_mean_oe_2(5));
    M_TDX_J2_mean_2   = TDX_J2_mean_oe_2(6);
    ex_TDX_J2_mean_2  = e_TDX_J2_mean_2*cos(omega_TDX_J2_mean_2);
    ey_TDX_J2_mean_2  = e_TDX_J2_mean_2*sin(omega_TDX_J2_mean_2);
    u_TDX_J2_mean_2   = M_TDX_J2_mean_2 + omega_TDX_J2_mean_2;
    TDX_J2_mean_2(idx,:)= [a_TDX_J2_mean_2, ex_TDX_J2_mean_2, ey_TDX_J2_mean_2, ...
                           i_TDX_J2_mean_2, RAAN_TDX_J2_mean_2, u_TDX_J2_mean_2];

    rel_J2_osc_qns_2(idx,:)= a_TSX_J2_2*compute_roes(a_TSX_J2_2, i_TSX_J2_2, e_TSX_J2_2, ...
                                          RAAN_TSX_J2_2, u_TSX_J2_2, ...
                                          a_TDX_J2_2, i_TDX_J2_2, e_TDX_J2_2, ...
                                          RAAN_TDX_J2_2, u_TDX_J2_2);
    rel_J2_mean_qns_2(idx,:)=a_TSX_J2_mean_2*compute_roes(a_TSX_J2_mean_2, i_TSX_J2_mean_2, e_TSX_J2_mean_2, ...
                                          RAAN_TSX_J2_mean_2, u_TSX_J2_mean_2, ...
                                          a_TDX_J2_mean_2, i_TDX_J2_mean_2, e_TDX_J2_mean_2, ...
                                          RAAN_TDX_J2_mean_2, u_TDX_J2_mean_2);

end

%–– labeling ––
elem_labels = { ...
  'a [m]',      ...
  'e_x [-]',        ...
  'e_y [-]',        ...
  'i [rad]',    ...
  'RAAN [rad]', ...
  'u [rad]'     ...
};
rel_labels = { ...
  'a\deltaa [m]',       ...
  'a\delta\lambda [m]', ...
  'a\deltae_x [m]',     ...
  'a\deltae_y [m]',     ...
  'a\deltai_x [m]',     ...
  'a\deltai_y [m]'      ...
};

case_names = { ...
  'Set 1 – no J2', ...
  'Set 1 – with J2', ...
  'Set 2 – no J2', ...
  'Set 2 – with J2'  ...
};

%–– collect your data in cell arrays ––
TSX_osc   = { TSX_osc_1,    TSX_J2_osc_1,    TSX_osc_2,    TSX_J2_osc_2    };
TSX_mean  = { TSX_mean_1,   TSX_J2_mean_1,   TSX_mean_2,   TSX_J2_mean_2   };
TDX_osc   = { TDX_osc_1,    TDX_J2_osc_1,    TDX_osc_2,    TDX_J2_osc_2    };
TDX_mean  = { TDX_mean_1,   TDX_J2_mean_1,   TDX_mean_2,   TDX_J2_mean_2   };
REL_osc   = { rel_osc_qns_1, rel_J2_osc_qns_1, rel_osc_qns_2, rel_J2_osc_qns_2 };
REL_mean  = { rel_mean_qns_1,rel_J2_mean_qns_1,rel_mean_qns_2,rel_J2_mean_qns_2 };

%–– main plotting loop ––
for k = 1:4
  % TSX absolute elements
  figure;
  for i = 1:6
    subplot(3,2,i);
    plot( t_orbit, TSX_osc{k}(:,i), '-', ...
          t_orbit, TSX_mean{k}(:,i), '--' );
    xlabel('Orbits'); ylabel(elem_labels{i});
    legend('Osculating','Mean','Location','best');
  end
  sgtitle([case_names{k} ' — TSX']);

  % TDX absolute elements
  figure;
  for i = 1:6
    subplot(3,2,i);
    plot( t_orbit, TDX_osc{k}(:,i), '-', ...
          t_orbit, TDX_mean{k}(:,i), '--' );
    xlabel('Orbits'); ylabel(elem_labels{i});
    legend('Osculating','Mean','Location','best');
  end
  sgtitle([case_names{k} ' — TDX']);

  % Relative QNS elements
  figure;
  for i = 1:6
    subplot(3,2,i);
    plot( t_orbit, REL_osc{k}(:,i), '-', ...
          t_orbit, REL_mean{k}(:,i), '--' );
    xlabel('Orbits'); ylabel(rel_labels{i});
    legend('Osculating','Mean','Location','best');
  end
  sgtitle([case_names{k} ' — Relative QNS']);
end

%Part d) Plot in RTN frame

TSX_noJ2 = { TSX_rv_1,       TSX_rv_2       };
TDX_noJ2 = { TDX_rv_1,       TDX_rv_2       };
TSX_J2   = { TSX_rv_J2_1,    TSX_rv_J2_2    };
TDX_J2   = { TDX_rv_J2_1,    TDX_rv_J2_2    };
set_names = {'Set 1', 'Set 2'};

for s = 1:2
  % preallocate
  rel_noJ2 = zeros(num_points,3);
  rel_J2   = zeros(num_points,3);

  % compute RTN-relative position for each case
  for idx = 1:num_points
    state1 = eci2rtn(TSX_noJ2{s}(idx,:)', TDX_noJ2{s}(idx,:)');
    rel_noJ2(idx,:) = state1(1:3)';   % [R, T, N]
    state2 = eci2rtn(TSX_J2{s}(idx,:)',   TDX_J2{s}(idx,:)');
    rel_J2(idx,:)   = state2(1:3)';
  end

  % 1) 2D projections
  figure;
  % TR plane: x = T, y = R
  subplot(1,3,1)
  plot(rel_noJ2(:,2), rel_noJ2(:,1), '-', ...
       rel_J2(:,2),   rel_J2(:,1),   '--')
  axis equal; grid on
  xlabel('T [m]'); ylabel('R [m]')
  title(['TR — ' set_names{s}])
  legend('no J2','J2','Location','best')

  % NR plane: x = N, y = R
  subplot(1,3,2)
  plot(rel_noJ2(:,3), rel_noJ2(:,1), '-', ...
       rel_J2(:,3),   rel_J2(:,1),   '--')
  axis equal; grid on
  xlabel('N [m]'); ylabel('R [m]')
  title(['NR — ' set_names{s}])

  % TN plane: x = T, y = N
  subplot(1,3,3)
  plot(rel_noJ2(:,2), rel_noJ2(:,3), '-', ...
       rel_J2(:,2),   rel_J2(:,3),   '--')
  axis equal; grid on
  xlabel('T [m]'); ylabel('N [m]')
  title(['TN — ' set_names{s}])

  % 2) 3D RTN scatter
  figure;
  plot3(rel_noJ2(:,1), rel_noJ2(:,2), rel_noJ2(:,3), '-', ...
        rel_J2(:,1),   rel_J2(:,2),   rel_J2(:,3),   '--')
  axis equal; grid on
  xlabel('R [m]'); ylabel('T [m]'); zlabel('N [m]')
  title(['3D RTN — ' set_names{s}])
  legend('no J2','J2','Location','best')
end

% Part e)

rel_osc_all  = { rel_osc_qns_1,  rel_J2_osc_qns_1,  rel_osc_qns_2,  rel_J2_osc_qns_2  };
rel_mean_all = { rel_mean_qns_1, rel_J2_mean_qns_1, rel_mean_qns_2, rel_J2_mean_qns_2 };
case_names   = { ...
  'Set 1 – no J2', ...
  'Set 1 – with J2', ...
  'Set 2 – no J2', ...
  'Set 2 – with J2'  ...
};

for k = 1:4
  osc  = rel_osc_all{k};
  mean = rel_mean_all{k};

  % 1) Relative eccentricity vector (dex vs dey)
  figure;
  plot( osc(:,3),  osc(:,4),  '-',  ...
        mean(:,3), mean(:,4), '--' );
  axis equal; grid on;
  xlabel('a\deltae_x [m]'); ylabel('a\deltae_y [m]');
  title([case_names{k} ' — Rel. Eccentricity Vector']);
  legend('Osculating','Mean','Location','best');

  % 2) Relative inclination vector (dix vs diy)
  figure;
  plot( osc(:,5),  osc(:,6),  '-',  ...
        mean(:,5), mean(:,6), '--' );
  axis equal; grid on;
  xlabel('a\deltai_x [m]'); ylabel('a\deltai_y [m]');
  title([case_names{k} ' — Rel. Inclination Vector']);
  legend('Osculating','Mean','Location','best');

  % 3) Relative mean longitude vs semi-major axis (dλ vs da)
  figure;
  plot( osc(:,2),  osc(:,1),  '-',  ...
        mean(:,2), mean(:,1), '--' );
  axis equal; grid on;
  xlabel('a\delta\lambda [m]'); ylabel('a\deltaa [m]');
  title([case_names{k} ' — d\lambda vs d a']);
  legend('Osculating','Mean','Location','best');
end

%Part f) Change on initial relative orbital elements & maneuver

%based on the dlambda' equation from Lecture 5/6 
% dlambda' = -21/2(gamma*sin(2i)*dix+1/7*da)
% for no drift in dlambda, we want dix = 0
% We change dix to 0 degrees
% To subtract 30 degrees from inclination, we use:
% dv = sqrt( 2v^2 - 4v*cos(i) )
% direction will be anti normal at u=0 or normal at u=180

a_TDX_f =  6886536.686;
e_TDX_f = 0.0001269;
i_TDX_f = deg2rad(97.4554);
RAAN_TDX_f = deg2rad(351.0106);
omega_TDX_f = deg2rad(-100.5043);
M_TDX_f = wrapTo2Pi(2*pi - omega_TDX_f);
u_TDX_f = omega_TDX_f + M_TDX_f %should be 0

TDX_oe_f = [a_TDX_f, e_TDX_f, i_TDX_f, RAAN_TDX_f, omega_TDX_f, M_TDX_f];
TDX_rv_f = oe2rv(TDX_oe_f, mu);

v = norm(TDX_rv_f(4:6));
dv = sqrt(2*v^2 - 4*v*cosd(30))


% Part g) Propagation with new qns roe
rel_qns_g       = [0, 100, 50, 100, 0, 200];
TDX_init_oe_g   = qns2oe(TSX_init_oe_2, rel_qns_g);
TDX_init_osc_oe_g = TDX_init_oe_g;
TDX_init_rv_g   = oe2rv(TDX_init_osc_oe_g, mu);

[t_out, TDX_rv_g]    = ode4(@compute_rates_rv_unperturbed, [tstart, tend]', TDX_init_rv_g, dt);
[t_out, TDX_rv_J2_g] = ode4(@compute_rates_rv_perturbed,   [tstart, tend]', TDX_init_rv_g, dt);

rtn_out_g      = zeros(num_points,3);
rtn_out_J2_g   = zeros(num_points,3);
TDX_osc_oe_g   = zeros(num_points,6);
TDX_mean_oe_g  = zeros(num_points,6);
rel_osc_g      = zeros(num_points,6);
rel_mean_g     = zeros(num_points,6);

for idx = 1:num_points

    state1 = eci2rtn(TSX_rv_2(idx,:)', TDX_rv_g(idx,:)');
    rtn_out_g(idx,:) = state1(1:3)';
    state2 = eci2rtn(TSX_rv_J2_2(idx,:)', TDX_rv_J2_g(idx,:)');
    rtn_out_J2_g(idx,:) = state2(1:3)';

    TSX_params_J2_2   = rv2oe(TSX_rv_J2_2(idx,:), mu);
    TDX_params_J2_2   = rv2oe(TDX_rv_J2_g(idx,:), mu);

    a_TSX_J2_2       = TSX_params_J2_2(1);
    e_TSX_J2_2       = TSX_params_J2_2(2);
    i_TSX_J2_2       = TSX_params_J2_2(3);
    RAAN_TSX_J2_2    = TSX_params_J2_2(4);
    omega_TSX_J2_2   = TSX_params_J2_2(5);
    nu_TSX_J2_2      = TSX_params_J2_2(6);
    e_vec_TSX_J2_2   = TSX_params_J2_2(11:12);
    M_TSX_J2_2       = true2mean(nu_TSX_J2_2, e_TSX_J2_2);
    u_TSX_J2_2       = M_TSX_J2_2 + omega_TSX_J2_2;
    TSX_J2_osc_2(idx,:)= [a_TSX_J2_2, e_vec_TSX_J2_2(:)', i_TSX_J2_2, ...
                          RAAN_TSX_J2_2, u_TSX_J2_2];

    TSX_J2_mean_oe_2  = osc2mean([a_TSX_J2_2, e_TSX_J2_2, i_TSX_J2_2, ...
                                 RAAN_TSX_J2_2, omega_TSX_J2_2, M_TSX_J2_2],1)'; %1 for J2
    a_TSX_J2_mean_2   = TSX_J2_mean_oe_2(1);
    e_TSX_J2_mean_2   = TSX_J2_mean_oe_2(2);
    i_TSX_J2_mean_2   = TSX_J2_mean_oe_2(3);
    RAAN_TSX_J2_mean_2= wrapTo2Pi(TSX_J2_mean_oe_2(4));
    omega_TSX_J2_mean_2=wrapTo2Pi(TSX_J2_mean_oe_2(5));
    M_TSX_J2_mean_2   = TSX_J2_mean_oe_2(6);
    ex_TSX_J2_mean_2  = e_TSX_J2_mean_2*cos(omega_TSX_J2_mean_2);
    ey_TSX_J2_mean_2  = e_TSX_J2_mean_2*sin(omega_TSX_J2_mean_2);
    u_TSX_J2_mean_2   = M_TSX_J2_mean_2 + omega_TSX_J2_mean_2;
    TSX_J2_mean_2(idx,:)= [a_TSX_J2_mean_2, ex_TSX_J2_mean_2, ey_TSX_J2_mean_2, ...
                           i_TSX_J2_mean_2, RAAN_TSX_J2_mean_2, u_TSX_J2_mean_2];

    a_TDX_J2_2       = TDX_params_J2_2(1);
    e_TDX_J2_2       = TDX_params_J2_2(2);
    i_TDX_J2_2       = TDX_params_J2_2(3);
    RAAN_TDX_J2_2    = TDX_params_J2_2(4);
    omega_TDX_J2_2   = TDX_params_J2_2(5);
    nu_TDX_J2_2      = TDX_params_J2_2(6);
    e_vec_TDX_J2_2   = TDX_params_J2_2(11:12);
    M_TDX_J2_2       = true2mean(nu_TDX_J2_2, e_TDX_J2_2);
    u_TDX_J2_2       = M_TDX_J2_2 + omega_TDX_J2_2;
    TDX_osc_oe_g(idx,:)= [a_TDX_J2_2, e_vec_TDX_J2_2(:)', i_TDX_J2_2, ...
                          RAAN_TDX_J2_2, u_TDX_J2_2];

    TDX_J2_mean_oe_2  = osc2mean([a_TDX_J2_2, e_TDX_J2_2, i_TDX_J2_2, ...
                                 RAAN_TDX_J2_2, omega_TDX_J2_2, M_TDX_J2_2],1)'; %1 for J2
    a_TDX_J2_mean_2   = TDX_J2_mean_oe_2(1);
    e_TDX_J2_mean_2   = TDX_J2_mean_oe_2(2);
    i_TDX_J2_mean_2   = TDX_J2_mean_oe_2(3);
    RAAN_TDX_J2_mean_2= wrapTo2Pi(TDX_J2_mean_oe_2(4));
    omega_TDX_J2_mean_2=wrapTo2Pi(TDX_J2_mean_oe_2(5));
    M_TDX_J2_mean_2   = TDX_J2_mean_oe_2(6);
    ex_TDX_J2_mean_2  = e_TDX_J2_mean_2*cos(omega_TDX_J2_mean_2);
    ey_TDX_J2_mean_2  = e_TDX_J2_mean_2*sin(omega_TDX_J2_mean_2);
    u_TDX_J2_mean_2   = M_TDX_J2_mean_2 + omega_TDX_J2_mean_2;
    TDX_mean_oe_g(idx,:)= [a_TDX_J2_mean_2, ex_TDX_J2_mean_2, ey_TDX_J2_mean_2, ...
                           i_TDX_J2_mean_2, RAAN_TDX_J2_mean_2, u_TDX_J2_mean_2];

    rel_osc_g(idx,:)= a_TSX_J2_2*compute_roes(a_TSX_J2_2, i_TSX_J2_2, e_TSX_J2_2, ...
                                          RAAN_TSX_J2_2, u_TSX_J2_2, ...
                                          a_TDX_J2_2, i_TDX_J2_2, e_TDX_J2_2, ...
                                          RAAN_TDX_J2_2, u_TDX_J2_2);
    rel_mean_g(idx,:)=a_TSX_J2_mean_2*compute_roes(a_TSX_J2_mean_2, i_TSX_J2_mean_2, e_TSX_J2_mean_2, ...
                                          RAAN_TSX_J2_mean_2, u_TSX_J2_mean_2, ...
                                          a_TDX_J2_mean_2, i_TDX_J2_mean_2, e_TDX_J2_mean_2, ...
                                          RAAN_TDX_J2_mean_2, u_TDX_J2_mean_2);

end

%–– Pull out each relative QNS component as before
ada_osc_g         = rel_osc_g(:,1);
adlambda_osc_g    = rel_osc_g(:,2);
adex_osc_g        = rel_osc_g(:,3);
adey_osc_g        = rel_osc_g(:,4);
adix_osc_g        = rel_osc_g(:,5);
adiy_osc_g        = rel_osc_g(:,6);

ada_mean_g        = rel_mean_g(:,1);
adlambda_mean_g   = rel_mean_g(:,2);
adex_mean_g       = rel_mean_g(:,3);
adey_mean_g       = rel_mean_g(:,4);
adix_mean_g       = rel_mean_g(:,5);
adiy_mean_g       = rel_mean_g(:,6);

figure;
plot(adex_osc_g, adey_osc_g, '-'   , 'DisplayName','Osculating');
hold on;
plot(adex_mean_g, adey_mean_g, '--', 'DisplayName','Mean');
axis equal;
xlabel('a\deltae_x [m]');
ylabel('a\deltae_y [m]');
legend('Location','best');
grid on;

figure;
plot(adix_osc_g, adiy_osc_g, '-'   , 'DisplayName','Osculating');
hold on;
plot(adix_mean_g, adiy_mean_g, '--', 'DisplayName','Mean');
axis equal;
xlabel('a\deltai_x [m]');
ylabel('a\deltai_y [m]');
legend('Location','best');
grid on;

figure;
plot(adlambda_osc_g, ada_osc_g,'-'   , 'DisplayName','Osculating');
hold on;
plot(adlambda_mean_g, ada_mean_g,'--', 'DisplayName','Mean');
axis equal;
ylabel('a\deltaa [m]');
xlabel('a\delta\lambda [m]');
legend('Location','best');
grid on;

R_g   = rtn_out_g(:,1);
T_g   = rtn_out_g(:,2);
N_g   = rtn_out_g(:,3);

R_J2 = rtn_out_J2_g(:,1);
T_J2 = rtn_out_J2_g(:,2);
N_J2 = rtn_out_J2_g(:,3);

figure;

% TR plane
subplot(1,3,1);
plot(T_g,   R_g,   '-',  T_J2,   R_J2,   '--');
xlabel('T [m]');    ylabel('R [m]');
axis equal;         grid on;
title('TR Plane');
legend('Unperturbed','J2 Perturbed');

% NR plane
subplot(1,3,2);
plot(N_g,   R_g,   '-',  N_J2,   R_J2,   '--');
xlabel('N [m]');    ylabel('R [m]');
axis equal;         grid on;
title('NR Plane');
legend('Unperturbed','J2 Perturbed');

% TN plane
subplot(1,3,3);
plot(T_g,   N_g,   '-',  T_J2,   N_J2,   '--');
xlabel('T [m]');    ylabel('N [m]');
axis equal;         grid on;
title('TN Plane');
legend('Unperturbed','J2 Perturbed');

figure;
plot3( R_g,  T_g,  N_g,  '-',  R_J2,  T_J2,  N_J2,  '--');
xlabel('R [m]');    ylabel('T [m]');    zlabel('N [m]');
axis equal;         grid on;
title('3D RTN Trajectories');
legend('Unperturbed','J2 Perturbed');


%Part h)

roes_2_h = zeros(num_points,6);
roes_6_h = zeros(num_points,6);
roes_2_h(1,:) = rel_qns;
roes_6_h(1,:) = rel_qns_g;

for idx = 2:num_points

    oe = TSX_J2_mean_2(idx,:);
    phi = stm_qns_j2(dt, oe);
    roes_2_h(idx,:) = (phi*roes_2_h(idx-1,:)');
    roes_6_h(idx,:) = (phi*roes_6_h(idx-1,:)');

end

figure;
plot(roes_2_h(:,3),  roes_2_h(:,4),          '-'); hold on;
plot(rel_J2_osc_qns_2(:,3), rel_J2_osc_qns_2(:,4), '--');
plot(rel_J2_mean_qns_2(:,3),rel_J2_mean_qns_2(:,4), '-.');
axis equal; grid on;
xlabel('a\delta e_x [m]');
ylabel('a\delta e_y [m]');
legend('STM‐Propagated ROE (2)', 'J_2 Osculating QNS (6)', 'J_2 Mean QNS (6)','Location','best');

figure;
plot(roes_6_h(:,3),  roes_6_h(:,4),          '-'); hold on;
plot(rel_osc_g(:,3), rel_osc_g(:,4),         '--');
plot(rel_mean_g(:,3),rel_mean_g(:,4),        '-.');
axis equal; grid on;
xlabel('a\delta e_x [m]');
ylabel('a\delta e_y [m]');
legend('STM‐Propagated ROE (6)', 'Osculating QNS (6)', 'Mean QNS (6)','Location','best');

figure;
plot(roes_2_h(:,5),  roes_2_h(:,6),          '-'); hold on;
plot(rel_J2_osc_qns_2(:,5), rel_J2_osc_qns_2(:,6), '--');
plot(rel_J2_mean_qns_2(:,5),rel_J2_mean_qns_2(:,6), '-.');
axis equal; grid on;
xlabel('a\delta i_x [m]');
ylabel('a\delta i_y [m]');
legend('STM‐Propagated ROE (2)', 'Osculating QNS (2)', 'J_2 Mean QNS (2)','Location','best');

figure;
plot(roes_6_h(:,5),  roes_6_h(:,6),          '-'); hold on;
plot(rel_osc_g(:,5), rel_osc_g(:,6),         '--');
plot(rel_mean_g(:,5),rel_mean_g(:,6),        '-.');
axis equal; grid on;
xlabel('a\delta i_x [m]');
ylabel('a\delta i_y [m]');
legend('STM‐Propagated ROE (6)', 'Osculating QNS (6)', 'Mean QNS (6)','Location','best');

figure;
plot(roes_2_h(:,2),  roes_2_h(:,1),          '-'); hold on;
plot(rel_J2_osc_qns_2(:,2), rel_J2_osc_qns_2(:,1), '--');
plot(rel_J2_mean_qns_2(:,2),rel_J2_mean_qns_2(:,1), '-.');
axis equal; grid on;
xlabel('a\delta\lambda [m]');
ylabel('a\delta a [m]');
legend('STM‐Propagated ROE (2)', ' Osculating QNS (2)', 'Mean QNS (2)','Location','best');

figure;
plot(roes_6_h(:,2),  roes_6_h(:,1),          '-'); hold on;
plot(rel_osc_g(:,2), rel_osc_g(:,1),         '--');
plot(rel_mean_g(:,2),rel_mean_g(:,1),        '-.');
axis equal; grid on;
xlabel('a\delta\lambda [m]');
ylabel('a\delta a [m]');
legend('STM‐Propagated ROE (6)', 'Osculating QNS (6)', 'Mean QNS (6)','Location','best');