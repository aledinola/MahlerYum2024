function [sol,sol_cpu,mom] = solve_model(Params,Gr,flag,vfoptions,simoptions,FnsToEval)

% Unpack variables from struct "flag":
do_cpu=flag.do_cpu;

% Unpack variables from struct "Gr":
% -- Dimensions
n_d = Gr.n_d;
n_a = Gr.n_a;
n_semiz = Gr.n_semiz;
n_z = Gr.n_z;
n_e = Gr.n_e;
N_j = Gr.N_j;
N_i = Gr.N_i;
% -- Grids
d_grid     = Gr.d_grid;
a_grid     = Gr.a_grid;
z_grid     = Gr.z_grid;
e_grid     = Gr.e_grid;
semiz_grid = vfoptions.semiz_grid;
theta_grid = Params.theta_i;
% -- Probability distributions
pi_z  = Gr.pi_z;
pi_e  = Gr.pi_e;
dist_health = Gr.dist_health;
% -- Toolkit options
PTypeDistParamNames = Gr.PTypeDistParamNames;

DiscParNam={'beta','sj'};

ReturnFn = @(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss,sigma,delta,nu,gamma,iota_j,psi)...
    Mod_ReturnFn(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss,sigma,delta,nu,gamma,iota_j,psi);

%% VFI: given parameters and grids, compute V and Policy
disp('ValueFnIter on GPU')

tic;
[V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,...
    z_grid,pi_z,ReturnFn,Params,DiscParNam,vfoptions);
t_vfi_toolkit=toc;

disp('Check size of value function')
disp('n_a,   n_semiz, n_z, n_e, Age')
disp(size(V.low))
disp([n_a,n_semiz, n_z,vfoptions.n_e,N_j])

disp('Check size of policy function')
% 3 = policy for d1(labor supply),d2(effort),a'(savings)
disp('   3,   n_a,  n_semiz, n_z, n_e, Age')
disp(size(Policy.low))
disp([3,n_a,n_semiz, n_z,vfoptions.n_e,N_j])

% VFI on the cpu, my code
if do_cpu==1
    disp('ValueFnIter on CPU')
    tic
    [V_cpu,Policy_cpu] = fun_vfi_cpu(Params,Gr,vfoptions);
    t_vfi_cpu=toc;

    % Convert Robert's notation to my notation
    [V_gpu,~,Policy_gpu] = reshape_VandPolicy(V,[],Policy,n_a,n_semiz,n_z,n_e,N_i,N_j);

    errV = max(abs(V_gpu-V_cpu),[],"all");
    errP = max(abs(Policy_gpu-Policy_cpu),[],"all");
    fprintf('errV: %f \n',errV)
    fprintf('errP: %f \n',errP)
end

%% Distribution

%% Initial distribution of agents at birth (j=1)
% Before we plot the life-cycle profiles we have to define how agents are
% at age j=1. We will give them all zero assets.
jequaloneDist=zeros([n_a,n_semiz,n_z,n_e],'gpuArray'); % Put no households anywhere on grid
% At j=1, agents have zero assets, non-degenerate distrib of zh, median zn,
% median medical spending shock, distrib of permanent types
jequaloneDist(1,:,floor((n_z+1)/2),floor((n_e+1)/2)) = dist_health;

AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age

%% We now compute the 'stationary distribution' of households
% default is parallel over e, 1 will loop over e, 2 will parfor loop over e
%simoptions.loopovere = 0;
tic
StatDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
t_dist_toolkit=toc;

disp('Check size of distribution')
disp('n_a,   n_semiz, n_z, n_e, Age')
disp(size(StatDist.low))
disp([n_a,n_semiz,n_z,vfoptions.n_e,N_j])
disp(size(StatDist.high))
disp([n_a,n_semiz,n_z,vfoptions.n_e,N_j])

if do_cpu==1
    tic
    %fun_StatDist_cpu(Params,Gr,Policy,joneDist,simoptions)
    StatDist_cpu = fun_StatDist_cpu(Params,Gr,Policy,jequaloneDist,simoptions);
    t_dist_cpu=toc;

    % Convert Robert's notation to my notation
    [~,StatDist_gpu,~] = reshape_VandPolicy(V,StatDist,Policy,n_a,n_semiz,n_z,n_e,N_i,N_j);

    err_D = max(abs(StatDist_gpu-StatDist_cpu),[],"all");
    fprintf('err_D: %f \n',err_D)
end

test_top = zeros(numel(N_i),1);
mass_PT = zeros(n_a,numel(N_i));
for ii=1:numel(N_i)
    myname = N_i{ii};
    mass_PT(:,ii) = sum(StatDist.(myname),[2,3,4,5]); %(a,zh,zn,e,agej)
    test_top(ii,1) = sum(mass_PT(end-10:end,ii));
    if max(test_top)>1e-4
        warning('Too much mass at top of the asset grid!') 
    end
end

%% Running times
disp('RUNNING TIMES')
fprintf('VFI toolkit       = %f \n',t_vfi_toolkit)
fprintf('StatDist toolkit  = %f \n',t_dist_toolkit)
if do_cpu==1
    fprintf('VFI cpu           = %f \n',t_vfi_cpu)
    fprintf('StatDist cpu      = %f \n',t_dist_cpu)
end
%fprintf('Moments (my code) = %f \n',time_mom)

%% Pack results into structures
sol = pack_into_struct(V,Policy,StatDist);

if do_cpu==1
    sol_cpu = pack_into_struct(V_cpu,Policy_cpu,StatDist_cpu);
else
    sol_cpu = [];
end

%% Moments

if nargout>2
    % Calculate model moments, conditional on age, health, etc.
    tic
    mom = fun_model_moments(Params,FnsToEval,V,StatDist,Policy,Gr,flag,simoptions);
    time_mom=toc;
    fprintf('Moments (my code) = %f \n',time_mom)

end

end %end function