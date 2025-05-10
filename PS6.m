clc; clear; close all;
addpath('mean_osc');
format long g;

global tol Re J2 mu s_d
tol = 10e-10;
Re = 6378137; % m
J2 = 1.082626e-3;
mu = 3.986004418e14; % (m^3/s^2)
s_d = 86400; % seconds per day

scenario = input( ...
    ['Select reconfiguration scenario:\n' ...
     '  1: M-D2 to M-D3\n' ...
     '  2: M-D3 to M-D4\n' ...
     'Your choice: '] );

switch scenario
    case 1
        % M-D2 to M-D3
        rel_qns_pre = [0, 0, 0, 300, 0, 400];
        rel_qns_post  = [0, 0, 0, 300, 0, 500];
        scenario_name = 'M-D2 to M-D3';
        
    case 2
        % M-D3 to M-D4
        rel_qns_pre = [0, 0, 0, 300, 0, 500];
        rel_qns_post  = [0, 0, 0, 500, 0, 300];
        scenario_name = 'M-D3 to M-D4';
        
    otherwise
        error('Invalid scenario');
end

fprintf('Running scenario "%s": rel_qns_init = [%g %g %g %g %g %g]\n rel_qns_fin = [%g %g %g %g %g %g]\n', ...
        scenario_name, [rel_qns_pre, rel_qns_post]);

% initial elements
a_TSX_init    = 6886536.686;  % m
e_TSX_init    = 0.0001264;
i_TSX_init    = deg2rad(97.4453);
RAAN_TSX_init = deg2rad(351.0108);
omega_TSX_init= deg2rad(101.2452);
M_TSX_init    = deg2rad(11.6520);
nu_TSX_init   = mean2true(M_TSX_init, e_TSX_init, tol);
u_TSX_init    = M_TSX_init + omega_TSX_init;

TSX_init_oe = [a_TSX_init, e_TSX_init, i_TSX_init, ...
               RAAN_TSX_init, omega_TSX_init, M_TSX_init];
TSX_init_rv = oe2rv(TSX_init_oe, mu);

[t_out, TSX_oe] = ode4(@compute_rates_GVE_J2, [tstart, tend]', TSX_init_oe', dt);
TSX_oe(:,6) = wrapTo2Pi(TSX_oe(:,6));

TDX_init_oe = qns2oe(TSX_init_oe, rel_qns_pre);
TDX_init_rv = oe2rv(TDX_init_oe, mu);
TDX_init_rtn = eci2rtn(TSX_init_rv, TDX_init_rv);
rel_state_init = [TDX_init_rtn; TSX_init_rv];

% timing parameters
tstart     = 0.0;
n          = sqrt(mu/a_TSX_init^3);
T          = 2*pi/n;
n_orbit    = 15;
tend       = n_orbit*T;
num_points = 1000;
dt         = (tend - tstart)/(num_points-1);
t_grid     = linspace(tstart, tend, num_points).';
t_orbit    = t_grid / T;
orbit_num = floor(t_orbit) + 1;

%delta V calculations


%propagation with delta Vs
state_out = zeros(num_points, 12);
state_out(1,:) = rel_state_init';

TDX_oe = zeros(num_points, 6);
TDX_oe(1,:) = TDX_init_oe;

rel_oe = zeros(num_points, 6);
rel_oe(1,:) = compute_roes(TSX_init_oe, TDX_init_oe);

state_cur = rel_state_init;

for idx = 2:num_points

    t_cur = t_grid(idx-1);
    t_next = t_grid(idx);

    TSX_ECI_cur = state_cur(7:12);
    TDX_RTN_cur = state_cur(1:6);

    % apply ΔV’s at the target indices


    % write the modified RTN back into state_cur
    state_cur(1:6)   = TDX_RTN_cur;

    [TDX_ECI_cur, ~] = rtn2eci(TSX_ECI_cur, TDX_RTN_cur);

    [~, state_next] = ode4(@compute_rates_rv_rel_unperturbed_RTN, [t_cur, t_next]', state_cur,  dt);
    state_next = state_next(2,:)';
    state_out(idx,:) = state_next;

    TSX_ECI_next = state_next(7:12);
    TDX_RTN_next = state_next(1:6);
    [TDX_ECI_next, R] = rtn2eci(TSX_ECI_next, TDX_RTN_next);
    
    a1 = TSX_oe(idx,1); e1 = TSX_oe(idx,2); i1 = TSX_oe(idx,3);
    RAAN1 = TSX_oe(idx,4); omega1 = TSX_oe(idx,5);
    M1 = wrapTo2Pi(TSX_oe(idx,6)); u1 = M1 + omega1;

    TDX_params = rv2oe(TDX_ECI_next, mu);
    TDX_oe(idx,:) = [TDX_params(1:5), true2mean(TDX_params(6), TDX_params(2))];
    
    a2 = TDX_oe(idx,1); e2 = TDX_oe(idx,2); i2 = TDX_oe(idx,3);
    RAAN2 = TDX_oe(idx,4); omega2 = TDX_oe(idx,5);
    M2 = wrapTo2Pi(TDX_oe(idx,6)); u2 = M2 + omega2;
    
    rel_oe(idx,:) = a1 * compute_roes([a1, e1, i1, RAAN1, omega1, M1], ...
                                      [a2, e2, i2, RAAN2, omega2, M2])';

    state_cur = state_next;

end