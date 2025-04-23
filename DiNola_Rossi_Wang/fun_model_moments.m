function [mom]=fun_model_moments(Params,FnsToEval,V_in,StatDist_in,Policy_in,Gr,flag,simoptions)
% OUTPUTS
% Age = j, model period
%   ave:     Average values for each variable, unconditional, scalar
%   ave_age: Average values for each variable, conditional on j, (1,J)
%   ave_health: Average values, conditional health (n_zh,1)
%   ave_age_health: Average values, conditional on j and on health (n_zh,J)

% Unpack params
Jr = Params.Jr;
% Unpack grids
n_a = Gr.n_a;
n_semiz = Gr.n_semiz;
n_z = Gr.n_z;
n_d = Gr.n_d;
n_e = Gr.n_e;
N_i = Gr.N_i;
N_j = Gr.N_j;
d_grid = Gr.d_grid;
a_grid = Gr.a_grid;
e_grid = Gr.e_grid;
z_grid = gather(Gr.z_grid); %dim: (n_z,1)
semiz_grid = Gr.zh_grid;
if iscell(N_i)
    n_theta = numel(N_i);
else
    n_theta = N_i;
end

if ~isequal(size(z_grid),[n_z,1])
    error('fun_model_moments: Size of z_grid.low is not correct')
end
if ~isequal(size(semiz_grid),[n_semiz,1])
    error('fun_model_moments: Size of z_grid.low is not correct')
end

%% mu is (a,semiz,z,e,theta,age)
%% Pol is (3,a,semiz,z,e,theta,age), integer-valued d1, d2 and a'
[V,mu,~] = reshape_VandPolicy(V_in,StatDist_in,Policy_in,n_a,n_semiz,n_z,n_e,N_i,N_j);

if any(mu<0,"all") || any(mu>1,"all")
    error('Distrib. mu must have all elements in [0,1]')
end
if abs(sum(mu,"all")-1)>1e-8
    error('Distrib. mu must sum to 1')
end

% Marginal distribution over age
mu_age = reshape(sum(mu,[1,2,3,4,5]),[1,N_j]);

%% Values on grid
Values_toolkit=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(StatDist_in,Policy_in,FnsToEval,Params,n_d,n_a,n_z,N_j,N_i,d_grid, a_grid, z_grid, simoptions);
% See line 242
% See also http://discourse.vfitoolkit.com/t/values-on-grid-with-semi-exogenous-shocks/345
%Values_struct=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(StatDist_in,Policy_in,FnsToEval,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

Values = fun_values_on_grid(Params,Policy_in,FnsToEval,d_grid,a_grid,semiz_grid,z_grid,e_grid,n_d,n_a,n_semiz,n_z,n_e,N_i,N_j);

names = fieldnames(Values);
names_young = {'hours','effort','health','labincome'};

%% Convert Values_struct into struct of arrays of size (a,zh,zn,e,theta,age)
Values_toolkit_arr = struct();
for ii=1:numel(names)
    name_var = names{ii};
    Values_toolkit_arr.(name_var) = zeros(n_a,n_semiz,n_z,n_e,n_theta,N_j);
    for j_c = 1:N_j
        % low type
        Values_toolkit_arr.(name_var)(:,:,:,:,1,j_c) = Values_toolkit.(name_var).low(:,:,:,:,j_c);
        % high type
        Values_toolkit_arr.(name_var)(:,:,:,:,2,j_c) = Values_toolkit.(name_var).high(:,:,:,:,j_c);
    end
end

if ~isequal(size(Values.assets),size(mu))
    error('fun_model_moments: Size of Values.X and Distribution are not compatible')
end

% Check that Values and Values_toolkit_arr are identical

%% Means
%% Unconditional average
ave   = struct();
for ii=1:numel(names)
    name_var = names{ii};
    ave.(name_var) = sum(Values.(name_var).*mu,"all");
end

%% Unconditional average, non-retirees only
% mu is (a,zn,zh,e,theta,age)
mass_young = sum(mu(:,:,:,:,:,1:Jr-1),'all');
ave_young   = struct();
for ii=1:numel(names_young)
    name_var = names_young{ii};
    ave_young.(name_var) = sum(Values.(name_var)(:,:,:,:,:,1:Jr-1).*mu(:,:,:,:,:,1:Jr-1),"all")/mass_young;
end

%% ave_age: Average values for each variable, conditional on j, (1,J)
ave_age = struct();
for ii=1:numel(names)
    name_var = names{ii};
    ave_age.(name_var) = squeeze(sum(Values.(name_var).*mu,[1,2,3,4,5]))./squeeze(sum(mu,[1,2,3,4,5]));
    ave_age.(name_var) = reshape(ave_age.(name_var),[1,N_j]);
end

%%   ave_health: Average values, conditional health (n_semiz,1)
% dim: (a,zh,zn,e,theta,age)
ave_health = struct(); 
% Each field has dim: (n_zh,1)
for ii=1:numel(names)
    name_var = names{ii};
    ave_health.(name_var) = squeeze(sum(Values.(name_var).*mu,[1,3,4,5,6]))./squeeze(sum(mu,[1,3,4,5,6]));
    ave_health.(name_var) = reshape(ave_health.(name_var),[n_semiz,1]);
end

%% Average values, conditional health (n_zh,1), non-retirees only
% mu is (a,zn,zh,e,theta,age)
mass_young = sum(mu(:,:,:,:,:,1:Jr-1),[1,3,4,5,6]);
ave_health_young   = struct();
for ii=1:numel(names_young)
    name_var = names_young{ii};
    ave_health_young.(name_var) = squeeze(sum(Values.(name_var)(:,:,:,:,:,1:Jr-1).*mu(:,:,:,:,:,1:Jr-1),[1,3,4,5,6]))...
        ./squeeze(mass_young);
end


%%   ave_age_health: Average values, conditional on j and on health (n_semiz,J)
% dim: (a,zh,zn,e,theta,age)
ave_age_health = struct();
% Each field has dim: (n_semiz,J)
for ii=1:numel(names)
    name_var = names{ii};
    ave_age_health.(name_var) = squeeze(sum(Values.(name_var).*mu,[1,3,4,5]))./squeeze(sum(mu,[1,3,4,5]));
%    ave_age_health.(name_var) = reshape(ave_age_health.(name_var),[n_zh,1]);
end

%% Ginis
%% Unconditional Gini coefficient
% gini = fun_gini(pop,val,makeplot)
gini = struct();
mu_vec=reshape(mu,[n_a*n_semiz*n_z*n_e*n_theta*N_j,1]);
for ii=1:numel(names)
    name_var = names{ii};
    values_vec = reshape(Values.(name_var),[n_a*n_semiz*n_z*n_e*n_theta*N_j,1]);
    gini.(name_var) = fun_gini(mu_vec,values_vec);
end

%% Gini coefficient by age
% gini = fun_gini(pop,val,makeplot)
gini_age = struct();
for ii=1:numel(names)
    name_var = names{ii};
    gini_age.(name_var) = zeros(1,N_j);
    for jj = 1:N_j
        mu_vec=reshape(mu(:,:,:,:,:,jj),[n_a*n_semiz*n_z*n_e*n_theta,1]); 
        % Normalize to sum to 1
        mu_vec = mu_vec/sum(mu_vec);
        values_vec = reshape(Values.(name_var)(:,:,:,:,:,jj),[n_a*n_semiz*n_z*n_e*n_theta,1]);
        gini_age.(name_var)(jj) = fun_gini(mu_vec,values_vec);
    end
end

%% Gini coefficient by health
% gini = fun_gini(pop,val,makeplot)
gini_health = struct();
for ii=1:numel(names)
    name_var = names{ii};
    gini_health.(name_var) = zeros(n_semiz,1);
    for zh_c = 1:n_semiz
        mu_vec=reshape(mu(:,zh_c,:,:,:,:),[n_a*n_z*n_e*n_theta*N_j,1]); 
        % Normalize to sum to 1
        mu_vec = mu_vec/sum(mu_vec);
        values_vec = reshape(Values.(name_var)(:,zh_c,:,:,:,:),[n_a*n_z*n_e*n_theta*N_j,1]);
        gini_health.(name_var)(zh_c) = fun_gini(mu_vec,values_vec);
    end
end

%% Correlation b/w labincome and health status zh (non-retirees)
% Vectorize labincome, health and mu
values_vec_labincome = reshape(Values.labincome(:,:,:,:,:,1:Jr-1),[n_a*n_semiz*n_z*n_e*n_theta*(Jr-1),1]);
values_vec_health = reshape(Values.health(:,:,:,:,:,1:Jr-1),[n_a*n_semiz*n_z*n_e*n_theta*(Jr-1),1]);
mu_vec = reshape(mu(:,:,:,:,:,1:Jr-1),[n_a*n_semiz*n_z*n_e*n_theta*(Jr-1),1]);
mu_vec = mu_vec/sum(mu_vec);

% means are ave.labincome and ave.health.
% Compute stdev
StdDev_labincome=sqrt(sum(mu_vec.*((values_vec_labincome-ave_young.labincome).^2)));
StdDev_health=sqrt(sum(mu_vec.*((values_vec_health-ave_young.health).^2)));

Numerator=sum((values_vec_labincome-ave_young.labincome).*(values_vec_health-ave_young.health).*mu_vec);
corr.labincome_health=Numerator/(StdDev_labincome*StdDev_health);

%% Tax moments

% Average tax rate of the average earner
% ave.taxable_income is the mean taxable income
tax_mean = fun_tax(ave.taxable_income,Params.tau0,Params.tau1);
tax.atr_mean = tax_mean/ave.taxable_income;

% Government net revenue
tax.govnet = sum((Values.taxrev-Values.pen-Values.b_tr-Values.med_sub).*mu,"all");

%% Welfare
% Integral of value function "V" with distribution "mu"

nan_ind = isnan(V);
welfare = sum(V(~nan_ind).*mu(~nan_ind),"all");

if anynan(V)
    nan_mass = sum(mu(nan_ind));
    warning('Value function V contains NaN elements')
    warning('NaN elements have total mass of %f \n',nan_mass)
end


% Pack results into a structure
ave.welfare = welfare;
mom = pack_into_struct(ave,ave_young,ave_age,ave_health,ave_health_young,...
    ave_age_health,gini,gini_age,gini_health,tax,corr);

end %end function fun_model_moments
