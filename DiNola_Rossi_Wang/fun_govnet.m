function [govnet] = fun_govnet(Params,Gr,flag,vfoptions,simoptions,FnsToEval)
% fun_govnet computes the government budget residual, defined as
%  Residual = income taxes + social security taxes - pensions - welfare
%  transfers - medical transfers

n_a = Gr.n_a;
n_z = Gr.n_z;
n_semiz = Gr.n_semiz;
n_d = Gr.n_d;
n_e = Gr.n_e;
N_i = Gr.N_i;
N_j = Gr.N_j;
d_grid = Gr.d_grid;
a_grid = Gr.a_grid;
semiz_grid = Gr.zh_grid;
e_grid = Gr.e_grid;
z_grid = Gr.z_grid; %dim: (n_z,1)

% Step 1: Call solve_model to get distribution and policy functions
sol = solve_model(Params,Gr,flag,vfoptions,simoptions,FnsToEval);
V       = sol.V;
Policy  = sol.Policy;
StatDist= sol.StatDist;

[~,mu_cpu,~] = reshape_VandPolicy(V,StatDist,Policy,n_a,n_semiz,n_z,n_e,N_i,N_j);

% Step 2: Calculate net government tax revenues, G in the paper. G is
% defined as total tax revenues-pensions-welfare-medical_susbidy
names = {'taxrev','pen','b_tr','med_sub'};
n_names = numel(names);
for ii = 1:n_names
    FnsToEval1.(names{ii}) = FnsToEval.(names{ii});
end

% Inputs are GPU-version
%Values_struct=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(StatDist,Policy,...
%    FnsToEval1,Params,n_d,n_a,n_z,N_j,N_i,d_grid,a_grid,z_grid,simoptions);

Values = fun_values_on_grid(Params,Policy,FnsToEval1,d_grid,a_grid,semiz_grid,z_grid,e_grid,n_d,n_a,n_semiz,n_z,n_e,N_i,N_j);

% %% Convert Values_struct into struct of arrays of size (a,zn,zh,e,theta,age)
% Values = struct();
% for ii=1:n_names
%     name_var = names{ii};
%     Values.(name_var) = zeros(n_a,n_semiz,n_z,n_e,n_theta,N_j);
%     for j_c = 1:N_j
%         % low type
%         Values.(name_var)(:,:,:,:,1,j_c) = Values_struct.(name_var).low(:,:,:,:,j_c);
%         % high type
%         Values.(name_var)(:,:,:,:,2,j_c) = Values_struct.(name_var).high(:,:,:,:,j_c);
%     end
% end

govnet = sum((Values.taxrev-Values.pen-Values.b_tr-Values.med_sub).*mu_cpu,"all");


end %end function "fun_govnet"