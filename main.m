%% Life-Cycle Model with health shocks
% TODO: 
% - Survival prob must depend on health state
% - Adjustment cost of effort with shock
% - Tax function and social transfers
clear,clc,close all,format long g
% Laptop
%myf = 'C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab';
% Desktop
myf = 'C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab';
addpath(genpath(myf))
% Add cpu-specific functions
%addpath('cpu')

%% How does VFI Toolkit think about this?
% Ordering of variables: (n,f,aprime,a,h,z)
% 1 Endogenous state variable:           a, assets 
% 2 Stochastic exogenous state variable: z, income shock and health shock
% 1 Permanent type (low vs high):        theta
%   Age:                                 j

%% Set parameters
% Age 25 to 100
start_age = 25;
last_age  = 100;
Params.J  = (last_age-start_age+1)/2; % =38, Number of period in life-cycle

% Grid sizes to use
n_d     = [3,30]; % Decision variables d1,d2: Labor choice (0,PT,FT), health effort
n_a     = 101;    % endogenous state: Asset holdings 
n_semiz = 2;      % Semi-exogenous state: Health status (0=unhealthy,1=healthy,2=dead)
n_z     = 5;      % Exogenous state: labor productivity shock
n_e     = 10;     % iid shocks: Adjustment cost, labor productivity
N_j     = Params.J; % Number of periods in finite horizon

% Demographics
Params.agej = 1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.Jr   = 21;
work_age = (1:1:Params.Jr-1)';

%% Permanent types
% e, e,            educ:     low, high
% beta,         patience: low, high
% productivity, theta:    low, high
% eta,          health:   low, high
% 4 permanent types, each of them has 2 points --> in total, 2^4=16 types

Names_i = {'e0b0p0h0','e1b0p0h0','e0b1p0h0','e1b1p0h0',...
    'e0b0p1h0','e1b0p1h0','e0b1p1h0','e1b1p1h0',...
    'e0b0p0h1','e1b0p0h1','e0b1p0h1','e1b1p0h1',...
    'e0b0p1h1','e1b0p1h1','e0b1p1h1','e1b1p1h1'};

if numel(Names_i)~=16
    error('Number of PT must be equal to 16')
end

% Disutility of working shifter. It depends on age j and on health status (h=0,1)

% Healthy, h=1 
nu_h1_1  = 2.634;
nu_h1_8  = 1.666;    
nu_h1_13 = 1.278;
nu_h1_20 = 1.714;

Params.nu_h1 = spline([1;8;13;20],[nu_h1_1;nu_h1_8;nu_h1_13;nu_h1_20],work_age);
Params.nu_h1 = [Params.nu_h1;zeros(Params.J-Params.Jr+1,1)];

% Unhealthy, h=0
nu_h0_1  = 2.412;
nu_h0_8  = 1.813;    
nu_h0_13 = 1.391;
nu_h0_20 = 2.415;

Params.nu_h0 = spline([1;8;13;20],[nu_h0_1;nu_h0_8;nu_h0_13;nu_h0_20],work_age);
Params.nu_h0 = [Params.nu_h0;zeros(Params.J-Params.Jr+1,1)];

figure
plot(work_age,Params.nu_h0(1:Params.Jr-1))
hold on
plot(work_age,Params.nu_h1(1:Params.Jr-1))
legend('Unhealthy','Healthy')
xlabel('Working age, j=1,..,Jr')
sgtitle('Health-dependent shifter for disutility of work')

Params.nu_e = 0.807; % Work disutility for non-college educated

% --- Wage profile
% Non-college (e=0)
Params.zeta0_e0 = 0.899;
Params.zeta1_e0 = 0.0616;
Params.zeta2_e0 = -0.0025;

Params.zeta0_e1 = 1.165;
Params.zeta1_e1 = 0.0874;
Params.zeta2_e1 = -0.0029;

% Wage penalty for unhealthy (h=0)
Params.wp_e0 = 0.178; % non-college 
Params.wp_e1 = 0.145; % college

% Permanent type theta
Params.theta_low  = -0.29;
Params.theta_high = 0.29;
Params.theta_prob = [0.5,0.5];

% Now create age-dependent component of wage, which depends on age j and on educ (e=0,1)
% Note: lambda in the paper includes also the wage penalty if unhealthy
Params.lambda_e0 = Params.zeta0_e0*exp(Params.zeta1_e0*(Params.agej-1)+Params.zeta2_e0*(Params.agej-1).^2);
Params.lambda_e1 = Params.zeta0_e1*exp(Params.zeta1_e1*(Params.agej-1)+Params.zeta2_e1*(Params.agej-1).^2);

figure
plot(1:1:Params.Jr-1,Params.lambda_e0(1:Params.Jr-1))
hold on
plot(1:1:Params.Jr-1,Params.lambda_e1(1:Params.Jr-1))
legend('Non-college','College')
xlabel('Age')
title('Deterministic component of log(wages)')

% --- Health effort parameters

% Disutility of effort (age-dependent) for healthy (h=1) and for e=0 (non-college)
iota_h1_e0_1  = 0.146;
iota_h1_e0_12 = 0.560;    
iota_h1_e0_20 = 1.048;
iota_h1_e0_31 = 1.603;

Params.iota_h1_e0 = spline([1;12;20;31],[iota_h1_e0_1;iota_h1_e0_12;iota_h1_e0_20;iota_h1_e0_31],Params.agej');
Params.iota_h1_e0 = [Params.iota_h1_e0(1:31);Params.iota_h1_e0(31)*ones(Params.J-31,1)];

% Disutility of effort (age-dependent) for unhealthy (h=0) and for e=0 (non-college)
iota_h0_e0_1  = 0.628;
iota_h0_e0_12 = 1.366;    
iota_h0_e0_20 = 1.650;
iota_h0_e0_31 = 0.735;

Params.iota_h0_e0 = spline([1;12;20;31],[iota_h0_e0_1;iota_h0_e0_12;iota_h0_e0_20;iota_h0_e0_31],Params.agej');
Params.iota_h0_e0 = [Params.iota_h0_e0(1:31);Params.iota_h0_e0(31)*ones(Params.J-31,1)];

% Disutility of effort (age-dependent) for healthy (h=1) for e=1 (college)
iota_h1_e1_1  = 0.0913;
iota_h1_e1_12 = 0.302;    
iota_h1_e1_20 = 0.740;
iota_h1_e1_31 = 1.366;

Params.iota_h1_e1 = spline([1;12;20;31],[iota_h1_e1_1;iota_h1_e1_12;iota_h1_e1_20;iota_h1_e1_31],Params.agej');
Params.iota_h1_e1 = [Params.iota_h1_e1(1:31);Params.iota_h1_e1(31)*ones(Params.J-31,1)];

% Disutility of effort (age-dependent) for unhealthy (h=0) for e=1 (college)
iota_h0_e1_1  = 0.469;
iota_h0_e1_12 = 0.997;    
iota_h0_e1_20 = 1.654;
iota_h0_e1_31 = 1.089;

Params.iota_h0_e1 = spline([1;12;20;31],[iota_h0_e1_1;iota_h0_e1_12;iota_h0_e1_20;iota_h0_e1_31],Params.agej');
Params.iota_h0_e1 = [Params.iota_h0_e1(1:31);Params.iota_h0_e1(31)*ones(Params.J-31,1)];

figure
subplot(1,2,1)
plot(Params.agej,Params.iota_h0_e0)
hold on
plot(Params.agej,Params.iota_h1_e0)
legend('Unhealthy','Healthy')
xlabel('Age, j')
title('Non-college')
subplot(1,2,2)
plot(Params.agej,Params.iota_h0_e1)
hold on
plot(Params.agej,Params.iota_h1_e1)
legend('Unhealthy','Healthy')
xlabel('Age, j')
title('College')
sgtitle('Effort disutility shifter, \iota')

% Adjustment cost for effort
Params.csi0 = 0.00012; 
Params.csi1 = 0.145;
Params.Bj = ones(1,N_j);
for jj=1:N_j
    Params.Bj(jj) = Params.csi0*exp(Params.csi1*(Params.agej(jj)-1));
end

% --- Remaining parameters
Params.r     = 1.04^2-1; % Net interest rate (one period = 2 years)
Params.kappa = 0.872;  % Shifter in utility of consumption (if h=0)
Params.sigma = 2;      % CRRA in utility of consumption
Params.b     = 13.11;  % Utility constant
Params.rho   = 0.975;  % Labor productivity autocorrelation
Params.sigx  = 0.0289; % Labor productivity standard deviation
Params.omega = 0.359;  % Pension scale
Params.gamma = 1;      % Elasticity of labor supply
Params.psi   = 1.115;  % Elasticity in cost of health effort function

% Tax function parameters
Params.tau_s = 0.321;
Params.tau_p = 0.128;

% --- Survival probabilities: depend on age j, on fixed education type
% e=0,1 and on health state h=0,1
sj_e1 = importdata("surv_CL.txt"); % e=1, college
sj_e0 = importdata("surv_HS.txt"); % e=0, non-college

Params.sj_e0 = sj_e0(:,1);
Params.sj_e1 = sj_e1(:,1);

figure
plot(Params.agej,Params.sj_e0)
hold on
plot(Params.agej,Params.sj_e1)
legend('Non-college','College')
xlabel('Age, j')
title('Conditional survival probabilities')

% --- Discount factor
mu_beta = 0.943;
delta_beta = 0.0284;
Params.beta_low  = mu_beta-delta_beta;
Params.beta_high = mu_beta+delta_beta;

% --- Health technology (probability of being in good health)
Params.pi0     = -0.905;
Params.lambda1 = 0.693;
Params.delta   = 2.311;
Params.gamma1  = 0.238;
Params.gamma2  = 0.632;

%% Set grids

% --- Assets grid
a_min = 0;
a_max = 30;
a_grid  = a_min+(a_max-a_min)*(linspace(0,1,n_a).^3)';

% --- Labor supply and health effort grids (d variables)
d1_grid = [0,0.5,1]'; % Labor supply
f_min = 0; % Minimal possible health effort
f_max = 1; % Maximum possible health effort
d2_grid = linspace(f_min,f_max,n_d(2))'; % Health effort
d_grid  = [d1_grid;d2_grid];

% --- Health status, semi-exogenous shock
semiz_grid = [0,1]'; % 0=unhealthy, 1=healthy

% --- Markov shock z to labor productivity
[z_grid,pi_z]=discretizeAR1_Rouwenhorst(0.0,Params.rho,Params.sigx,n_z);

% Labor earnings in last period before retirement. It is the same for
% everyone and depends only on education fixed type. Used to determine
% pension benefits
z_med = z_grid(floor((n_z+1)/2));
Params.y_pen_e0 = d1_grid(3)*Params.theta_prob(1)*exp(Params.lambda_e0(Params.Jr-1)+Params.theta_low+z_med)+...
    Params.theta_prob(2)*exp(Params.lambda_e0(Params.Jr-1)+Params.theta_high+z_med);
Params.y_pen_e1 = d1_grid(3)*Params.theta_prob(1)*exp(Params.lambda_e1(Params.Jr-1)+Params.theta_low+z_med)+...
    Params.theta_prob(2)*exp(Params.lambda_e1(Params.Jr-1)+Params.theta_high+z_med);

% --- IID shocks: adjustment cost and transitory shock eps to labor prod
% Adjustment cost is uniform random variable in [0,B(j)]
e_grid = zeros(n_e,N_j);
pi_e   = zeros(n_e,N_j);
for jj=1:N_j
    e_grid(:,jj) = linspace(0,Params.Bj(jj),n_e)';
    pi_e(:,jj)   = 1/n_e(1); %each point has same prob
end

%% Set all parameters that depend on permanent type

% Names_i = {'e0b0p0h0','e1b0p0h0','e0b1p0h0','e1b1p0h0',...
%     'e0b0p1h0','e1b0p1h0','e0b1p1h0','e1b1p1h0',...
%     'e0b0p0h1','e1b0p0h1','e0b1p0h1','e1b1p0h1',...
%     'e0b0p1h1','e1b0p1h1','e0b1p1h1','e1b1p1h1'};

% Cost of health effort, iota, depends on age j, on fixed educ type e=0,1
% and on health state h=0,1
% -- here for h=0 (unhealthy)
Params.iota_h0.e0b0p0h0 = Params.iota_h0_e0;
Params.iota_h0.e1b0p0h0 = Params.iota_h0_e1;
Params.iota_h0.e0b1p0h0 = Params.iota_h0_e0;
Params.iota_h0.e1b1p0h0 = Params.iota_h0_e1;
Params.iota_h0.e0b0p1h0 = Params.iota_h0_e0;
Params.iota_h0.e1b0p1h0 = Params.iota_h0_e1;
Params.iota_h0.e0b1p1h0 = Params.iota_h0_e0;
Params.iota_h0.e1b1p1h0 = Params.iota_h0_e1;
Params.iota_h0.e0b0p0h1 = Params.iota_h0_e0;
Params.iota_h0.e1b0p0h1 = Params.iota_h0_e1;
Params.iota_h0.e0b1p0h1 = Params.iota_h0_e0;
Params.iota_h0.e1b1p0h1 = Params.iota_h0_e1;
Params.iota_h0.e0b0p1h1 = Params.iota_h0_e0;
Params.iota_h0.e1b0p1h1 = Params.iota_h0_e1;
Params.iota_h0.e0b1p1h1 = Params.iota_h0_e0;
Params.iota_h0.e1b1p1h1 = Params.iota_h0_e1;
% -- here for h=1 (healthy)
Params.iota_h1.e0b0p0h0 = Params.iota_h1_e0;
Params.iota_h1.e1b0p0h0 = Params.iota_h1_e1;
Params.iota_h1.e0b1p0h0 = Params.iota_h1_e0;
Params.iota_h1.e1b1p0h0 = Params.iota_h1_e1;
Params.iota_h1.e0b0p1h0 = Params.iota_h1_e0;
Params.iota_h1.e1b0p1h0 = Params.iota_h1_e1;
Params.iota_h1.e0b1p1h0 = Params.iota_h1_e0;
Params.iota_h1.e1b1p1h0 = Params.iota_h1_e1;
Params.iota_h1.e0b0p0h1 = Params.iota_h1_e0;
Params.iota_h1.e1b0p0h1 = Params.iota_h1_e1;
Params.iota_h1.e0b1p0h1 = Params.iota_h1_e0;
Params.iota_h1.e1b1p0h1 = Params.iota_h1_e1;
Params.iota_h1.e0b0p1h1 = Params.iota_h1_e0;
Params.iota_h1.e1b0p1h1 = Params.iota_h1_e1;
Params.iota_h1.e0b1p1h1 = Params.iota_h1_e0;
Params.iota_h1.e1b1p1h1 = Params.iota_h1_e1;


% Deterministic component of wage, for e=0 and e=1, and for age j
Params.lambda.e0b0p0h0 = Params.lambda_e0;
Params.lambda.e1b0p0h0 = Params.lambda_e1;
Params.lambda.e0b1p0h0 = Params.lambda_e0;
Params.lambda.e1b1p0h0 = Params.lambda_e1;
Params.lambda.e0b0p1h0 = Params.lambda_e0;
Params.lambda.e1b0p1h0 = Params.lambda_e1;
Params.lambda.e0b1p1h0 = Params.lambda_e0;
Params.lambda.e1b1p1h0 = Params.lambda_e1;
Params.lambda.e0b0p0h1 = Params.lambda_e0;
Params.lambda.e1b0p0h1 = Params.lambda_e1;
Params.lambda.e0b1p0h1 = Params.lambda_e0;
Params.lambda.e1b1p0h1 = Params.lambda_e1;
Params.lambda.e0b0p1h1 = Params.lambda_e0;
Params.lambda.e1b0p1h1 = Params.lambda_e1;
Params.lambda.e0b1p1h1 = Params.lambda_e0;
Params.lambda.e1b1p1h1 = Params.lambda_e1;

% Permanent productivity type for wage, theta
Params.theta.e0b0p0h0 = Params.theta_low;
Params.theta.e1b0p0h0 = Params.theta_low;
Params.theta.e0b1p0h0 = Params.theta_low;
Params.theta.e1b1p0h0 = Params.theta_low;
Params.theta.e0b0p1h0 = Params.theta_high;
Params.theta.e1b0p1h0 = Params.theta_high;
Params.theta.e0b1p1h0 = Params.theta_high;
Params.theta.e1b1p1h0 = Params.theta_high;
Params.theta.e0b0p0h1 = Params.theta_low;
Params.theta.e1b0p0h1 = Params.theta_low;
Params.theta.e0b1p0h1 = Params.theta_low;
Params.theta.e1b1p0h1 = Params.theta_low;
Params.theta.e0b0p1h1 = Params.theta_high;
Params.theta.e1b0p1h1 = Params.theta_high;
Params.theta.e0b1p1h1 = Params.theta_high;
Params.theta.e1b1p1h1 = Params.theta_high;

% Survival probability, depends on age j and on e=0,1 (non-college,college)
Params.sj.e0b0p0h0 = Params.sj_e0;
Params.sj.e1b0p0h0 = Params.sj_e1;
Params.sj.e0b1p0h0 = Params.sj_e0;
Params.sj.e1b1p0h0 = Params.sj_e1;
Params.sj.e0b0p1h0 = Params.sj_e0;
Params.sj.e1b0p1h0 = Params.sj_e1;
Params.sj.e0b1p1h0 = Params.sj_e0;
Params.sj.e1b1p1h0 = Params.sj_e1;
Params.sj.e0b0p0h1 = Params.sj_e0;
Params.sj.e1b0p0h1 = Params.sj_e1;
Params.sj.e0b1p0h1 = Params.sj_e0;
Params.sj.e1b1p0h1 = Params.sj_e1;
Params.sj.e0b0p1h1 = Params.sj_e0;
Params.sj.e1b0p1h1 = Params.sj_e1;
Params.sj.e0b1p1h1 = Params.sj_e0;
Params.sj.e1b1p1h1 = Params.sj_e1;

% Discount factor, depends on permanent beta type
Params.beta.e0b0p0h0 = Params.beta_low;
Params.beta.e1b0p0h0 = Params.beta_low;
Params.beta.e0b1p0h0 = Params.beta_high;
Params.beta.e1b1p0h0 = Params.beta_high;
Params.beta.e0b0p1h0 = Params.beta_low;
Params.beta.e1b0p1h0 = Params.beta_low;
Params.beta.e0b1p1h0 = Params.beta_high;
Params.beta.e1b1p1h0 = Params.beta_high;
Params.beta.e0b0p0h1 = Params.beta_low;
Params.beta.e1b0p0h1 = Params.beta_low;
Params.beta.e0b1p0h1 = Params.beta_high;
Params.beta.e1b1p0h1 = Params.beta_high;
Params.beta.e0b0p1h1 = Params.beta_low;
Params.beta.e1b0p1h1 = Params.beta_low;
Params.beta.e0b1p1h1 = Params.beta_high;
Params.beta.e1b1p1h1 = Params.beta_high;

% Pension benefits, depends only on education fixed type e=0,1
% (non-college, college)
Params.y_pen.e0b0p0h0 = Params.y_pen_e0;
Params.y_pen.e1b0p0h0 = Params.y_pen_e1;
Params.y_pen.e0b1p0h0 = Params.y_pen_e0;
Params.y_pen.e1b1p0h0 = Params.y_pen_e1;
Params.y_pen.e0b0p1h0 = Params.y_pen_e0;
Params.y_pen.e1b0p1h0 = Params.y_pen_e1;
Params.y_pen.e0b1p1h0 = Params.y_pen_e0;
Params.y_pen.e1b1p1h0 = Params.y_pen_e1;
Params.y_pen.e0b0p0h1 = Params.y_pen_e0;
Params.y_pen.e1b0p0h1 = Params.y_pen_e1;
Params.y_pen.e0b1p0h1 = Params.y_pen_e0;
Params.y_pen.e1b1p0h1 = Params.y_pen_e1;
Params.y_pen.e0b0p1h1 = Params.y_pen_e0;
Params.y_pen.e1b0p1h1 = Params.y_pen_e1;
Params.y_pen.e0b1p1h1 = Params.y_pen_e0;
Params.y_pen.e1b1p1h1 = Params.y_pen_e1;

% Fixed education type e=0,1
Params.educ.e0b0p0h0 = 0;
Params.educ.e1b0p0h0 = 1;
Params.educ.e0b1p0h0 = 0;
Params.educ.e1b1p0h0 = 1;
Params.educ.e0b0p1h0 = 0;
Params.educ.e1b0p1h0 = 1;
Params.educ.e0b1p1h0 = 0;
Params.educ.e1b1p1h0 = 1;
Params.educ.e0b0p0h1 = 0;
Params.educ.e1b0p0h1 = 1;
Params.educ.e0b1p0h1 = 0;
Params.educ.e1b1p0h1 = 1;
Params.educ.e0b0p1h1 = 0;
Params.educ.e1b0p1h1 = 1;
Params.educ.e0b1p1h1 = 0;
Params.educ.e1b1p1h1 = 1;

% Fixed health type eta=0,1 (unhealthy, healthy)
Params.eta.e0b0p0h0 = 0;
Params.eta.e1b0p0h0 = 0;
Params.eta.e0b1p0h0 = 0;
Params.eta.e1b1p0h0 = 0;
Params.eta.e0b0p1h0 = 0;
Params.eta.e1b0p1h0 = 0;
Params.eta.e0b1p1h0 = 0;
Params.eta.e1b1p1h0 = 0;
Params.eta.e0b0p0h1 = 1;
Params.eta.e1b0p0h1 = 1;
Params.eta.e0b1p0h1 = 1;
Params.eta.e1b1p0h1 = 1;
Params.eta.e0b0p1h1 = 1;
Params.eta.e1b0p1h1 = 1;
Params.eta.e0b1p1h1 = 1;
Params.eta.e1b1p1h1 = 1;

% Wage penalty wp. Depends on education e=0,1 (non-college, college)
Params.wp.e0b0p0h0 = Params.wp_e0;
Params.wp.e1b0p0h0 = Params.wp_e1;
Params.wp.e0b1p0h0 = Params.wp_e0;
Params.wp.e1b1p0h0 = Params.wp_e1;
Params.wp.e0b0p1h0 = Params.wp_e0;
Params.wp.e1b0p1h0 = Params.wp_e1;
Params.wp.e0b1p1h0 = Params.wp_e0;
Params.wp.e1b1p1h0 = Params.wp_e1;
Params.wp.e0b0p0h1 = Params.wp_e0;
Params.wp.e1b0p0h1 = Params.wp_e1;
Params.wp.e0b1p0h1 = Params.wp_e0;
Params.wp.e1b1p0h1 = Params.wp_e1;
Params.wp.e0b0p1h1 = Params.wp_e0;
Params.wp.e1b0p1h1 = Params.wp_e1;
Params.wp.e0b1p1h1 = Params.wp_e0;
Params.wp.e1b1p1h1 = Params.wp_e1;


%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(n,f,aprime,a,h,z,agej,Jr,lambda,wp,theta,y_pen,r,kappa,nu_h0,nu_h1,iota_h0,iota_h1,gamma,psi,sigma,nu_e,b) ...
    MahlerYum2024_ReturnFn(n,f,aprime,a,h,z,agej,Jr,lambda,wp,theta,y_pen,r,kappa,nu_h0,nu_h1,iota_h0,iota_h1,gamma,psi,sigma,nu_e,b);

%% Set options for the toolkit
vfoptions.verbose    = 1; % Just so we can see feedback on progress

vfoptions.n_semiz    = n_semiz;
vfoptions.semiz_grid = semiz_grid;
% Define the transition probabilities of the semi-exogenous states
vfoptions.SemiExoStateFn=@(h,hprime,f,educ,eta,pi0,lambda1,delta,gamma1,gamma2) MahlerYum2024_SemiExoStateFn(h,hprime,f,educ,eta,pi0,lambda1,delta,gamma1,gamma2);
% It is hardcoded that only the 'last' decision variable can influence the transition probabilities of the semi-exogenous states
% The semi-exogenous states must be included in the return fn, fns to evaluate, etc. The ordering must be that semi-exogenous states come after this period endogenous state(s) and before any markov exogenous states, so (...,a,semiz,z,...)

% Now combine z and e together
% z_grid and pi_z will be passed directly as inputs
% e will be added to vfoptions and simoptions
vfoptions.n_e     = n_e;
vfoptions.e_grid  = e_grid;
vfoptions.pi_e    = pi_e;

% We also need to tell simoptions about the semi-exogenous states
simoptions.n_e            = vfoptions.n_e;
simoptions.e_grid         = vfoptions.e_grid;
simoptions.pi_e           = vfoptions.pi_e;
simoptions.n_semiz        = vfoptions.n_semiz;
simoptions.semiz_grid     = vfoptions.semiz_grid;
simoptions.SemiExoStateFn = vfoptions.SemiExoStateFn;

% %% Value function iteration
% disp('Test ValueFnIter')
% tic;
% [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
% toc
% 
% size(V.e0b0p0h0)
% disp([n_a,n_semiz,n_z,n_e,N_j])
% 
% size(Policy.e0b0p0h0)
% % which is the same as
% disp([length(n_d)+length(n_a),n_a,n_semiz,n_z,n_e,N_j])
% 
% % d1=labor, d2=health effort, a'
% pol_aprime_ind = squeeze(Policy.e1b1p1h1(3,:,:,:,:,:));
% pol_aprime = a_grid(pol_aprime_ind);
% max_aprime = max(pol_aprime,[],'all');
% 
% figure
% plot(a_grid,pol_aprime(:,:,1,1,1))

%% If you want to take a look at what the whole 'semi-exogenous transition matrix' looks like (it is created automatically by codes) it will look like
N_i = numel(Names_i);
pi_semiz_J=zeros([n_semiz,n_semiz,n_d(2),N_j]); % Note that endogenous state is the first, then the conditional transition matrix for shocks

for PT_c=1:N_i
for jj=1:N_j
for d2_c=1:n_d(2)
for hprime_c=1:n_semiz
for h_c=1:n_semiz
    % Note: the -1 turn the index into the value
    d2_val = d2_grid(d2_c);
    h_val = semiz_grid(h_c);
    hprime_val = semiz_grid(hprime_c);
    pi_semiz_J(h_c,hprime_c,d2_c,jj,PT_c)=MahlerYum2024_SemiExoStateFn(h_val,hprime_val,d2_val,Params.educ.(Names_i{PT_c}),Params.eta.(Names_i{PT_c}),Params.pi0,Params.lambda1,Params.delta,Params.gamma1,Params.gamma2);
end
end
end
end
end

% How does pi_semiz(2,1) = prob of getting sick changes with effort
ale = squeeze(pi_semiz_J(2,1,:,20,16));

% Make sure the 'rows' sum to one
for PT_c=1:N_i
for jj=1:N_j
for d2_c=1:n_d(2)
    temp=sum(pi_semiz_J(:,:,d2_c,jj,PT_c),2);
    if any(abs(temp-1)>10^(-14))
        temp-1
    end
end
end
end