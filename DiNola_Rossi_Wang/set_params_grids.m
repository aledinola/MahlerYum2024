function [Params,Gr,flag,vfoptions,simoptions,FnsToEval] = set_params_grids()

%% Set options
flag.do_plots = 0;
flag.results  = 'results';
flag.dir_name = 'model3'; % Name of subfolder in results (used in txt_export_all)
flag.do_bench = 1; % 1 = benchmark model, 0 otherwise
flag.do_cpu   = 0;

flag.saveG_fname = fullfile(flag.results,'govnet_bench.txt');

vfoptions.verbose = 1;
vfoptions.divideandconquer = 0;
vfoptions.level1n = 11;

simoptions.verbose = 0;
%simoptions.whichstats = [1,1,1,0,0,0,0];

if ~isfolder(flag.results)
    mkdir(flag.results)
end

%% Begin setting up to use VFI Toolkit to solve
% (d1,d2,aprime,a,semiz,z,e,..)
% Grid sizes to use
n_d     = [3,2]; % 1) d1:Labor supply, 2) d2:Health effort
n_a     = 101;   % Endogenous asset holdings
n_semiz = 2;     % Health status, semi-exog shock
n_z     = 7;     % Exogenous labor productivity shock
n_e     = 5;     % iid shock to medical expenses

%% Demographics

Params.agejshifter=20; % age=agej+agejshifter, where agej=1,..,J
Params.J=100-Params.agejshifter; % Total number of periods in life-cycle
% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.age=Params.agej+Params.agejshifter;
Params.Jr=45;
N_j = Params.J; % Number of periods in finite horizon

%% Permanent type PT
N_i = {'low','high'}; % No college vs college, or fixed income effect
n_theta = numel(N_i);
sigma_theta = 0.242; % Variance of PT
% Impose the two moment conditions:
% 0.5*theta(1)+0.5*theta(2)=0 (mean zero) ==> theta(2)=-theta(1)
% 0.5*theta(1)^2+0.5*theta(2)^2=sigma_theta  ==> theta(1)=sqrt(sigma_theta)
theta_grid = exp([-sqrt(sigma_theta),sqrt(sigma_theta)]);
Params.theta_i = theta_grid;
PTypeDistParamNames={'theta_dist'};
Params.theta_dist=[0.5,0.5]; 

% Discount rate
Params.beta = 0.98;
% Preferences
Params.sigma = 1.5; % Coeff of relative risk aversion in u(c)
Params.nu    = 1; % Weight on labor disutility 
Params.gamma = 1; % Elasticity of labor supply  
%Params.iota  = 1; % Weight on disutility health effort
Params.iota_j  = linspace(0,2,N_j);
Params.psi   = 1; % Elasticity health effort 

% Prices
Params.w=1;    % Wage
Params.r=0.04; % Interest rate

% Government. Net_income = tau0*gross_income^(1-tau1) 
Params.deduc = 0;    % if zero, can deduct 100% of medical expenses
Params.tau0  = 0.90; % Calibrate so that agg_taxes/agg_income=21%
Params.tau1  = 0.08; % Benchmark 0.08; flat 0; high prog. 0.16;
Params.sub   = 0.0;  % Subsidy for medical expenses (if =1, then no medical expenses)
Params.tau_ss = 0.12; % Social security tax rate
Params.cap_ss = 10000; % Social security cap

% Age-dependent labor productivity units
Params.kappa_j = zeros(1,Params.J);
Params.kappa_j(1:Params.Jr-1) = [1.0000, 1.0719, 1.1438, 1.2158, 1.2842, 1.3527, ...
              1.4212, 1.4897, 1.5582, 1.6267, 1.6952, 1.7217, ...
              1.7438, 1.7748, 1.8014, 1.8279, 1.8545, 1.8810, ...
              1.9075, 1.9341, 1.9606, 1.9623, 1.9640, 1.9658, ...
              1.9675, 1.9692, 1.9709, 1.9726, 1.9743, 1.9760, ...
              1.9777, 1.9700, 1.9623, 1.9546, 1.9469, 1.9392, ...
              1.9315, 1.9238, 1.9161, 1.9084, 1.9007, 1.8354, ...
              1.7701, 1.7048];

% Pensions
%pension_low=0.5*Params.w*Params.theta_i(1)*sum(Params.kappa_j)/(Params.Jr-1);
%pension_high=0.5*Params.w*Params.theta_i(2)*sum(Params.kappa_j)/(Params.Jr-1);
pension_low  = 0.5*Params.w*sum(Params.kappa_j)/(Params.Jr-1);
pension_high = pension_low;
Params.pen_j.low  = zeros(1,Params.J);
Params.pen_j.high = zeros(1,Params.J);
Params.pen_j.low(Params.Jr:Params.J) = pension_low;
Params.pen_j.high(Params.Jr:Params.J) = pension_high;

% Survival probability (must be of size [1,Params.J])
Params.sj = [1.00000, 0.99923, 0.99914, 0.99914, 0.99912, ...
                      0.99906, 0.99908, 0.99906, 0.99907, 0.99901, ...
                      0.99899, 0.99896, 0.99893, 0.99890, 0.99887, ...
                      0.99886, 0.99878, 0.99871, 0.99862, 0.99853, ...
                      0.99841, 0.99835, 0.99819, 0.99801, 0.99785, ...
                      0.99757, 0.99735, 0.99701, 0.99676, 0.99650, ...
                      0.99614, 0.99581, 0.99555, 0.99503, 0.99471, ...
                      0.99435, 0.99393, 0.99343, 0.99294, 0.99237, ...
                      0.99190, 0.99137, 0.99085, 0.99000, 0.98871, ...
                      0.98871, 0.98721, 0.98612, 0.98462, 0.98376, ...
                      0.98226, 0.98062, 0.97908, 0.97682, 0.97514, ...
                      0.97250, 0.96925, 0.96710, 0.96330, 0.95965, ...
                      0.95619, 0.95115, 0.94677, 0.93987, 0.93445, ...
                      0.92717, 0.91872, 0.91006, 0.90036, 0.88744, ...
                      0.87539, 0.85936, 0.84996, 0.82889, 0.81469, ...
                      0.79705, 0.78081, 0.76174, 0.74195, 0.72155, ...
                      0.00000];
Params.sj(end)=0; % The last period (j=J) value of sj is irrelevant

% Health effort: how it affects transitions
Params.pi0 = 1;
Params.pi1 = 0.5;

% Income penalty for bad health. It varies by permanent type. We assume
% that the low type loses more if bad health happens
Params.varrho.low  = exp(-0.3);
Params.varrho.high = exp(-0.3);

% Preference shifter for health: shifter(zh)=1-delta*zh
Params.delta = -1;%0.1;

% Consumption floor
Params.c_floor = 0.01;

%% Grids
a_min  = 0;
a_max  = 600;
a_spacing = 3;
a_grid = a_min+(a_max-a_min)*(linspace(0,1,n_a).^a_spacing)'; 

%% Decision variables d_grid
n_grid = [0,0.5,1.0]'; % Jung and Tran, JEEA
f_grid = [0,1]'; % no effort vs effort (for ENDO)
%f_grid = [0.302591,0.302591]'; % for EXOG
d_grid = [n_grid;f_grid];

%% z_grid and pi_z are exogenous shock, NOT age-dependent and NOT dependent on d
% z = zn
% e is iis shock to medical expenses, present only if zh=bad
Params.sigma_eps = sqrt(0.022); % Stdev of innovation to z, called eps
Params.rho       = 0.985;       % Persistence
Params.sig_e     = sqrt(0.2);   % Stdev of iid shock e
% z_grid_J: structure, each PT is size[n_z(1)+n_z(2),N_j]
% pi_z_J: structure, each PT is size[n_z(1)*n_z(2),N_j]
[p_bb_j,p_gb_j,zh_grid,dist_health,z_grid,pi_z,e_grid,pi_e]=create_shocks(n_semiz,n_z,n_e,N_j,Params);
% To use e variables we have to put them into the vfoptions and simoptions
vfoptions.n_e    =n_e;
vfoptions.e_grid =e_grid;
vfoptions.pi_e   =pi_e;
simoptions.n_e   =vfoptions.n_e;
simoptions.e_grid=vfoptions.e_grid;
simoptions.pi_e  =vfoptions.pi_e;

% Set up semi-exogenous shock parameters
Params.p_bb_j = p_bb_j;
Params.p_gb_j = p_gb_j;

% Set up the semi-exogneous states
vfoptions.n_semiz   =n_semiz;
vfoptions.semiz_grid=zh_grid;
% Define the transition probabilities of the semi-exogenous states
vfoptions.SemiExoStateFn=@(zh,zh_prime,d2_effort,p_bb_j,p_gb_j,pi0,pi1) Mod_SemiExoStateFn(zh,zh_prime,d2_effort,p_bb_j,p_gb_j,pi0,pi1);
% It is hardcoded that only the 'last' decision variable can influence the transition probabilities of the semi-exogenous states
% The semi-exogenous states must be included in the return fn, fns to evaluate, etc. The ordering must be that semi-exogenous states come after this period endogenous state(s) and before any markov exogenous states, so (...,a,semiz,z,...)

% We also need to tell simoptions about the semi-exogenous states
simoptions.n_semiz=vfoptions.n_semiz;
simoptions.semiz_grid=vfoptions.semiz_grid;
simoptions.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions.d_grid = d_grid;

% Out-of-pocket medical expenses, age-dependent. Actual medical expenses
% are defined as hc(j,e)=medspend_j(j)*e, where j=age and e is iid shock.
Params.medspend_j = linspace(0.3,0.9,Params.J);
%Params.medspend_j = zeros(1,Params.J);
hc = zeros(Params.J,length(e_grid));
for jj=1:Params.J
    hc(jj,:) = Params.medspend_j(jj)*e_grid;
end

% % Own calculation: simulate zh_prob for j=2,..J, starting from dist_health
% % zh_prob.PT has size (n_zh,age)
% zh_prob.low = zeros(n_z(2),Params.J);
% zh_prob.low(:,1) = dist_health;
% zh_prob.high = zeros(n_z(2),Params.J);
% zh_prob.high(:,1) = dist_health;
% for jj=1:Params.J-1
%     zh_prob.low(:,jj+1) = zh_prob.low(:,jj)'*pi_zh_J.low(:,:,jj);
%     zh_prob.high(:,jj+1) = zh_prob.high(:,jj)'*pi_zh_J.high(:,:,jj);
% end

% Marginal distribution of households over age
Params.mewj=ones(1,Params.J); 
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one

%% FnsToEvaluate and conditional restrictions

%--- Functions to be evaluatated for statistics, etc.
FnsToEval.cons = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss)...
    Mod_cons(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss);
FnsToEval.income = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss)...
    Mod_income(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,sub);
FnsToEval.labincome = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss)...
    Mod_labincome(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,varrho,w);
FnsToEval.assets = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) a; % a is the current asset holdings
FnsToEval.frac_badhealth=@(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) (zh==0); % indicator for z=0 (bad health) 
FnsToEval.taxrev = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss)...
    Mod_taxrev(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss);
FnsToEval.hours = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) n;
FnsToEval.effort = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) f;
FnsToEval.health = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) zh;
FnsToEval.medical = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) (zh==0)*medspend_j*e;
FnsToEval.med_sub = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) sub*(zh==0)*medspend_j*e;
FnsToEval.taxable_income= @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss)...
    Mod_taxable_income(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,sub);
FnsToEval.pen = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) pen_j;
FnsToEval.b_tr = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss)...
    Mod_b_tr(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,sub);

%--- Conditional restrictions. Must return either 0 or 1
condres.sick = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) (zh==0);
condres.healthy = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss) (zh==1);
% Add additional field to simoptions
simoptions.conditionalrestrictions = condres;

%% Check that Mod_SemiExoStateFn(zh,zh_prime,d2_effort,p_bb_j,p_gb_j,pi0,pi1)
% generates a plausible transition matrix
%               (semiz,semiz'|    f,theta,j)
pi_semiz_J=zeros([n_semiz,n_semiz,n_d(2),n_theta,N_j]); 
for j_c=1:N_j         % age
for theta_c=1:n_theta % PT
for d2_c=1:n_d(2)     % effort
    d2_val = f_grid(d2_c);
    for zh_c=1:n_semiz
        zh = zh_grid(zh_c);
        for zhprime_c=1:n_semiz
            zh_prime = zh_grid(zhprime_c);
            pi_semiz_J(zh_c,zhprime_c,d2_c,theta_c,j_c)=Mod_SemiExoStateFn(zh,zh_prime,d2_val,...
                Params.p_bb_j(j_c),Params.p_gb_j(theta_c,j_c),Params.pi0,Params.pi1);
        end
    end
end %end effort
end %end PT
end %end age
% Make sure the 'rows' sum to one
for j_c=1:N_j
    for theta_c=1:n_theta
        for d2_c=1:n_d(2)
            temp=sum(pi_semiz_J(:,:,d2_c,theta_c,j_c),2);
            if any(abs(temp-1)>10^(-14))
                warning('Row of pi_semiz_J does not sum to one!')
                fprintf('theta_c = %d, j_c = %d, d2_c = %d \n',theta_c,j_c,d2_c)
            end
        end
    end
end

%% Pack all grids and probabilities into a structure
Gr = pack_into_struct(n_d,d_grid,a_grid,n_a,n_semiz,n_z,n_e,N_i,N_j,...
    zh_grid,dist_health,z_grid,pi_z,e_grid,pi_e,PTypeDistParamNames,n_theta,...
    pi_semiz_J);

end %end function set_params_grids