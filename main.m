%% Life-Cycle Model with health shocks
% TODO: 
% - Survival prob must depend on health state
% - Adjustment cost of effort with shock
% - Tax function and social transfers
clear,clc,close all,format long g
% Laptop
myf = 'C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab';
% Desktop
%myf = 'C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab';
addpath(genpath(myf))
% Add cpu-specific functions
%addpath('cpu')

%% How does VFI Toolkit think about this?
% Ordering of variables: (n,f,aprime,a,h,z)
% 2 Decision variables: n (labor supply) and f (health effort)
% 1 Endogenous state variable:           aprime, a (assets) 
% 1 Semi-exogenous shock:                h (health status)
% 1 Exogenous shock:                     z (labor productivity)
% Age:                                   j
% 4 Permanent types:                    
%   educ  = fixed education
%   beta  = fixed patience
%   theta = fixed productivity
%   eta   = fixed health type

%% Set parameters

% --- Demographics
start_age = 25;
last_age  = 100;
Params.J  = (last_age-start_age+1)/2; % =38, Number of period in life-cycle
Params.agej = 1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.Jr   = 21;
work_age    = (1:1:Params.Jr-1)'; % 1,2,..,Jr-1 where Jr-1 = last working age

% --- Grid sizes to use
n_d     = [3,10];   % Decision variables d1,d2: Labor choice (0,PT,FT), health effort
n_a     = 301;      % endogenous state: Asset holdings 
n_semiz = 3;        % Semi-exogenous state: Health status (0=unhealthy,1=healthy,2=dead)
n_z     = 5;        % Exogenous state: labor productivity shock
N_j     = Params.J; % Number of periods in finite horizon

% --- Remaining parameters
Params.r     = 1.04^2-1; % Net interest rate (one period = 2 years)
Params.kappa = 0.872;    % Shifter in utility of consumption (if h=0)
Params.sigma = 2;        % CRRA in utility of consumption
Params.b     = 13.11;    % Utility constant
Params.rho   = 0.975;    % Labor productivity autocorrelation
Params.sigx  = 0.0289;   % Labor productivity standard deviation
Params.omega = 0.359;    % Pension scale
Params.gamma = 1;        % Elasticity of labor supply
Params.psi   = 1.115;    % Elasticity in cost of health effort function
Params.csi0  = 0.00012;  % Adjustment cost effort (not used)
Params.csi1  = 0.145;    % Adjustment cost effort (not used)

% Tax function parameters (HSV or Benabou functional form)
Params.tau_s = 0.321;    % Level parameter (roughly, average tax rate)
Params.tau_p = 0.128;    % Progressivity parameter

% Discount factor
mu_beta          = 0.943;
delta_beta       = 0.0284;
Params.beta_low  = mu_beta-delta_beta;
Params.beta_high = mu_beta+delta_beta;

% Part time and full time hours
Params.n_pt = 0.5;
Params.n_ft = 1.0;
 
% Health technology (probability of being in good health in j+1)
Params.pi0     = -0.905;
Params.lambda1 = 0.693;
Params.delta   = 2.311;
Params.gamma1  = 0.238;
Params.gamma2  = 0.632;

%% Permanent types
% e,            educ:     low, high
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
nu_h1_1  = 2.634;  % Age 1
nu_h1_8  = 1.666;  % Age 8   
nu_h1_13 = 1.278;  % Age 13
nu_h1_20 = 1.714;  % Age 20

Params.nu_h1 = spline([1;8;13;20],[nu_h1_1;nu_h1_8;nu_h1_13;nu_h1_20],work_age);
% Not relevant for age>=Jr, so we set arbitrary value of zero
Params.nu_h1 = [Params.nu_h1;zeros(Params.J-Params.Jr+1,1)];

% Unhealthy, h=0
nu_h0_1  = 2.412; % Age 1
nu_h0_8  = 1.813; % Age 8   
nu_h0_13 = 1.391; % Age 13
nu_h0_20 = 2.415; % Age 20

Params.nu_h0 = spline([1;8;13;20],[nu_h0_1;nu_h0_8;nu_h0_13;nu_h0_20],work_age);
Params.nu_h0 = [Params.nu_h0;zeros(Params.J-Params.Jr+1,1)];

figure
plot(work_age,Params.nu_h0(work_age))
hold on
plot(work_age,Params.nu_h1(work_age))
legend('Unhealthy','Healthy')
xlabel('Working age, j=1,..,Jr-1')
sgtitle('Health-dependent shifter for disutility of work')

Params.nu_e = 0.807; % Work disutility for non-college educated (e=0)

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

Params.iota_h1_e0 = spline([1;12;20;31],[iota_h1_e0_1;iota_h1_e0_12;iota_h1_e0_20;iota_h1_e0_31],(1:1:31)');
% Every age-dependent parameter must have J elements. I freeze at last one
Params.iota_h1_e0 = [Params.iota_h1_e0;Params.iota_h1_e0(31)*ones(Params.J-31,1)];

% Disutility of effort (age-dependent) for unhealthy (h=0) and for e=0 (non-college)
iota_h0_e0_1  = 0.628;
iota_h0_e0_12 = 1.366;    
iota_h0_e0_20 = 1.650;
iota_h0_e0_31 = 0.735;

Params.iota_h0_e0 = spline([1;12;20;31],[iota_h0_e0_1;iota_h0_e0_12;iota_h0_e0_20;iota_h0_e0_31],(1:1:31)');
Params.iota_h0_e0 = [Params.iota_h0_e0;Params.iota_h0_e0(31)*ones(Params.J-31,1)];

% Disutility of effort (age-dependent) for healthy (h=1) for e=1 (college)
iota_h1_e1_1  = 0.0913;
iota_h1_e1_12 = 0.302;    
iota_h1_e1_20 = 0.740;
iota_h1_e1_31 = 1.366;

Params.iota_h1_e1 = spline([1;12;20;31],[iota_h1_e1_1;iota_h1_e1_12;iota_h1_e1_20;iota_h1_e1_31],(1:1:31)');
Params.iota_h1_e1 = [Params.iota_h1_e1;Params.iota_h1_e1(31)*ones(Params.J-31,1)];

% Disutility of effort (age-dependent) for unhealthy (h=0) for e=1 (college)
iota_h0_e1_1  = 0.469;
iota_h0_e1_12 = 0.997;    
iota_h0_e1_20 = 1.654;
iota_h0_e1_31 = 1.089;

Params.iota_h0_e1 = spline([1;12;20;31],[iota_h0_e1_1;iota_h0_e1_12;iota_h0_e1_20;iota_h0_e1_31],(1:1:31)');
Params.iota_h0_e1 = [Params.iota_h0_e1;Params.iota_h0_e1(31)*ones(Params.J-31,1)];

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
xlabel('Age, j=1,..,J')
title('College')
sgtitle('Effort disutility shifter, \iota')

% --- Survival probabilities: depend on age j, on fixed education type
%     e=0,1 and on health state h=0,1
sj_e1 = importdata("surv_CL.txt"); % e=1 (college), col1 is h=1, col2 is h=0
sj_e0 = importdata("surv_HS.txt"); % e=0 (non-college), col1 is h=1, col2 is h=0 

Params.sj_h0_e0 = sj_e0(:,2); % unhealthy, non-college
Params.sj_h1_e0 = sj_e0(:,1); % healthy, non-college
Params.sj_h0_e1 = sj_e0(:,2); % unhealthy, college
Params.sj_h1_e1 = sj_e0(:,1); % healthy, college

figure
subplot(1,2,1)
plot(Params.agej,Params.sj_h0_e0)
hold on
plot(Params.agej,Params.sj_h1_e0,'--')
legend('Unhealthy','Healthy','Location','southwest')
xlabel('Age, j')
title('Non-college')

subplot(1,2,2)
plot(Params.agej,Params.sj_h0_e1)
hold on
plot(Params.agej,Params.sj_h1_e1,'--')
legend('Unhealthy','Healthy','Location','southwest')
xlabel('Age, j')
title('College')
sgtitle('Conditional survival probabilities')

%% Set grids

% --- Assets grid
a_min  = 0;
a_max  = 30;
a_grid = a_min+(a_max-a_min)*(linspace(0,1,n_a).^3)';

% --- Labor supply and health effort grids (d variables)
d1_grid = [0,Params.n_pt,Params.n_ft]'; % Labor supply
f_min   = 0; % Minimal possible health effort
f_max   = 1; % Maximum possible health effort
d2_grid = linspace(f_min,f_max,n_d(2))'; % Health effort
d_grid  = [d1_grid;d2_grid];

% --- Health status, semi-exogenous shock
semiz_grid = [0,1,2]'; % 0=unhealthy, 1=healthy, 2=dead

% --- Markov shock z to labor productivity
[z_grid,pi_z]=discretizeAR1_Rouwenhorst(0.0,Params.rho,Params.sigx,n_z);

% Labor earnings in last period before retirement. It is the same for
% everyone and depends only on education fixed type. Used to determine
% pension benefits
z_med = z_grid(floor((n_z+1)/2));
% e0 = non-college
Params.y_pen_e0 = Params.n_ft*Params.theta_prob(1)*exp(Params.lambda_e0(Params.Jr-1)+Params.theta_low+z_med)+...
    Params.theta_prob(2)*exp(Params.lambda_e0(Params.Jr-1)+Params.theta_high+z_med);
% e1s = college
Params.y_pen_e1 = Params.n_ft*Params.theta_prob(1)*exp(Params.lambda_e1(Params.Jr-1)+Params.theta_low+z_med)+...
    Params.theta_prob(2)*exp(Params.lambda_e1(Params.Jr-1)+Params.theta_high+z_med);

%% Set all parameters that depend on permanent type
% e=0,1 education
% b=0,1 patience
% p=0,1 fixed productivity (theta) 
% h=0,1 fixed health type (eta)
% Names_i = {'e0b0p0h0','e1b0p0h0','e0b1p0h0','e1b1p0h0',...
%     'e0b0p1h0','e1b0p1h0','e0b1p1h0','e1b1p1h0',...
%     'e0b0p0h1','e1b0p0h1','e0b1p0h1','e1b1p0h1',...
%     'e0b0p1h1','e1b0p1h1','e0b1p1h1','e1b1p1h1'};

% Cost of health effort, iota, depends on health state h=0,1 
% on fixed educ type e=0,1 and on age j 
% -- here for h=0 (unhealthy)
Params.iota_h0.e0b0p0h0 = Params.iota_h0_e0; % dim: (N_j,1)
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
Params.iota_h1.e0b0p0h0 = Params.iota_h1_e0; % dim: (N_j,1)
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

% Deterministic component of wage, depends on fixed educ type e=0 and e=1
% and on age j
Params.lambda.e0b0p0h0 = Params.lambda_e0; % dim: (N_j,1)
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

% Permanent productivity type (theta), depends on p=0,1
Params.theta.e0b0p0h0 = Params.theta_low; % dim: scalar
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

% Survival probability, depends on health state h=0,1 (not on fixed health type eta!) 
% on fixed educ type e=0,1 and on age j
% --- Here h0 (unhealthy)
Params.sj_h0.e0b0p0h0 = Params.sj_h0_e0; % dim: (N_j,1)
Params.sj_h0.e1b0p0h0 = Params.sj_h0_e1;
Params.sj_h0.e0b1p0h0 = Params.sj_h0_e0;
Params.sj_h0.e1b1p0h0 = Params.sj_h0_e1;
Params.sj_h0.e0b0p1h0 = Params.sj_h0_e0;
Params.sj_h0.e1b0p1h0 = Params.sj_h0_e1;
Params.sj_h0.e0b1p1h0 = Params.sj_h0_e0;
Params.sj_h0.e1b1p1h0 = Params.sj_h0_e1;
Params.sj_h0.e0b0p0h1 = Params.sj_h0_e0;
Params.sj_h0.e1b0p0h1 = Params.sj_h0_e1;
Params.sj_h0.e0b1p0h1 = Params.sj_h0_e0;
Params.sj_h0.e1b1p0h1 = Params.sj_h0_e1;
Params.sj_h0.e0b0p1h1 = Params.sj_h0_e0;
Params.sj_h0.e1b0p1h1 = Params.sj_h0_e1;
Params.sj_h0.e0b1p1h1 = Params.sj_h0_e0;
Params.sj_h0.e1b1p1h1 = Params.sj_h0_e1;

% --- Here h1 (healthy)
Params.sj_h1.e0b0p0h0 = Params.sj_h1_e0; % dim: (N_j,1)
Params.sj_h1.e1b0p0h0 = Params.sj_h1_e1;
Params.sj_h1.e0b1p0h0 = Params.sj_h1_e0;
Params.sj_h1.e1b1p0h0 = Params.sj_h1_e1;
Params.sj_h1.e0b0p1h0 = Params.sj_h1_e0;
Params.sj_h1.e1b0p1h0 = Params.sj_h1_e1;
Params.sj_h1.e0b1p1h0 = Params.sj_h1_e0;
Params.sj_h1.e1b1p1h0 = Params.sj_h1_e1;
Params.sj_h1.e0b0p0h1 = Params.sj_h1_e0;
Params.sj_h1.e1b0p0h1 = Params.sj_h1_e1;
Params.sj_h1.e0b1p0h1 = Params.sj_h1_e0;
Params.sj_h1.e1b1p0h1 = Params.sj_h1_e1;
Params.sj_h1.e0b0p1h1 = Params.sj_h1_e0;
Params.sj_h1.e1b0p1h1 = Params.sj_h1_e1;
Params.sj_h1.e0b1p1h1 = Params.sj_h1_e0;
Params.sj_h1.e1b1p1h1 = Params.sj_h1_e1;

% Discount factor, depends on permanent beta type b=0,1
Params.beta.e0b0p0h0 = Params.beta_low;   % dim: scalar
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
Params.wp.e0b0p0h0 = Params.wp_e0;  % dim: scalar
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
DiscountFactorParamNames={'beta'};

ReturnFn=@(n,f,aprime,a,h,z,agej,Jr,lambda,wp,theta,y_pen,r,kappa,nu_h0,nu_h1,iota_h0,iota_h1,gamma,psi,sigma,nu_e,b,educ) ...
    MahlerYum2024_ReturnFn(n,f,aprime,a,h,z,agej,Jr,lambda,wp,theta,y_pen,r,kappa,nu_h0,nu_h1,iota_h0,iota_h1,gamma,psi,sigma,nu_e,b,educ);

%% Set options for the toolkit
vfoptions.verbose    = 1; % Just so we can see feedback on progress

vfoptions.n_semiz    = n_semiz;
vfoptions.semiz_grid = semiz_grid;
% Define the transition probabilities of the semi-exogenous states
vfoptions.SemiExoStateFn=@(h,hprime,f,sj_h0,sj_h1,educ,eta,pi0,lambda1,delta,gamma1,gamma2) MahlerYum2024_SemiExoStateFn(h,hprime,f,sj_h0,sj_h1,educ,eta,pi0,lambda1,delta,gamma1,gamma2);
% It is hardcoded that only the 'last' decision variable can influence the transition probabilities of the semi-exogenous states
% The semi-exogenous states must be included in the return fn, fns to evaluate, etc. The ordering must be that semi-exogenous states come after this period endogenous state(s) and before any markov exogenous states, so (...,a,semiz,z,...)

% We also need to tell simoptions about the semi-exogenous states
simoptions.n_semiz        = vfoptions.n_semiz;
simoptions.semiz_grid     = vfoptions.semiz_grid;
simoptions.SemiExoStateFn = vfoptions.SemiExoStateFn;
simoptions.d_grid         = d_grid;

%% Value function iteration
disp('Test ValueFnIter')
tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
toc

size(V.e0b0p0h0)
disp([n_a,n_semiz,n_z,N_j])

size(Policy.e0b0p0h0)
% which is the same as
disp([length(n_d)+length(n_a),n_a,n_semiz,n_z,N_j])

% d1=labor, d2=health effort, a'
pol_aprime_ind = squeeze(Policy.e1b1p1h1(3,:,:,:,:));
pol_aprime = a_grid(pol_aprime_ind);
max_aprime = max(pol_aprime,[],'all');

figure
plot(a_grid,pol_aprime(:,:,1,1))

save temp.mat 

%% Initial distribution of agents at birth (j=1)

% We have to define how agents are at age j=1. We will give them all zero assets.
% State space: a,semiz,z,age,PT, where PT=educ,beta,theta,eta

% Distribution of permanent types (2^4 = 16 types in total)
PTypeDistParamNames={'PTypeDist'};

% Order of PTypes: eta,theta,beta,educ, i.e. eta varies first, etc.
PTypeDist_MY = [0.062, 0.133, 0.070, 0.101, 0.061, 0.099, 0.063, 0.102,...
    0.034, 0.033, 0.024, 0.045, 0.029, 0.047, 0.025, 0.072]';
PTypeDist = zeros(numel(PTypeDist_MY),1);

ind=1;
for eta_c=1:2
for theta_c=1:2
for beta_c=1:2
for e_c=1:2
    ind_MY = sub2ind([2,2,2,2],eta_c,theta_c,beta_c,e_c);
    PTypeDist(ind) = PTypeDist_MY(ind_MY);
    ind = ind+1;
end
end
end
end

Params.PTypeDist = PTypeDist;
jequaloneDist=zeros([n_a,n_semiz,n_z],'gpuArray'); 

% All agents start with zero assets, health status from data, and median productivity shock 
% At j=1, we asssume 10% unhealthy, 90% healthy and 0% dead
jequaloneDist(1,:,floor((n_z+1)/2)) = [0.1,0.9,0];

AgeWeightsParamNames={'ageweights'};
Params.ageweights=ones(N_j,1)/N_j;

StatDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

%% FnsToEvaluate
% n,f: Two decision variables d1 and d2
% aprime and a: Endogenous state (tomorrow and today)
% h: Semi-exogenous shock
% z: productivity shock
FnsToEvaluate.assets=@(n,f,aprime,a,h,z) a; % a is the current asset holdings

%--- Conditional restrictions. Must return either 0 or 1
condres.alive = @(n,f,aprime,a,h,z) (h==0 || h==1);
condres.dead  = @(n,f,aprime,a,h,z) (h==2);
% Add additional field to simoptions
simoptions.conditionalrestrictions = condres;

%% Calculate the life-cycle profiles

AgeStats=LifeCycleProfiles_FHorz_Case1_PType(StatDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
% By default, this includes both the 'grouped' statistics, like
% AgeConditionalStats.earnings.Mean
% Which are calculated across all permanent types of agents.
% And also the 'conditional on permanent type' statistics, like 
% AgeConditionalStats.earnings.patient.Mean

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
    hprime_val = semiz_grid(hprime_c);   %h,hprime,f,sj_h0,sj_h1,educ,eta,pi0,lambda1,delta,gamma1,gamma2
    pi_semiz_J(h_c,hprime_c,d2_c,jj,PT_c)=MahlerYum2024_SemiExoStateFn(h_val,hprime_val,d2_val,Params.sj_h0.(Names_i{PT_c})(jj),Params.sj_h1.(Names_i{PT_c})(jj),Params.educ.(Names_i{PT_c}),Params.eta.(Names_i{PT_c}),Params.pi0,Params.lambda1,Params.delta,Params.gamma1,Params.gamma2);
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
        warning('Rows do not sum to one!')
    end
    if ~isequal(pi_semiz_J(3,:,d2_c,jj,PT_c),[0,0,1])
        warning('Prob(h''|h=2) is not correct')
    end
end
end
end