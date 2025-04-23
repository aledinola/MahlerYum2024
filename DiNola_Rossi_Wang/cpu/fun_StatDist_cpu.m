function [StatDist] = fun_StatDist_cpu(Params,Gr,Policy_in,joneDist,simoptions)

joneDist = gather(joneDist);

% Unoack from structures
n_a     = Gr.n_a;
n_semiz = Gr.n_semiz;
n_z     = Gr.n_z;
n_e     = Gr.n_e;
N_i     = Gr.N_i;
N_j     = Gr.N_j;
pi_z    = gather(Gr.pi_z);
pi_e    = gather(Gr.pi_e);
pi_semiz_J = gather(Gr.pi_semiz_J);
mewj       = Params.mewj;
theta_dist = Params.theta_dist;

if iscell(N_i)
    N_i = numel(N_i);
end

% Convert policy from GPU to CPU and add dimension for PT
%  Policy has dim: (1:3,a,semiz,z,e,theta,agej)
[~,~,Policy] = reshape_VandPolicy([],[],Policy_in,n_a,n_semiz,n_z,n_e,N_i,N_j);

if ~isequal(size(Policy),[3,n_a,n_semiz,n_z,n_e,N_i,N_j])
    disp(size(Policy))
    error('fun_StatDist_cpu: Policy has wrong size')
end

%% Pre-allocate output array
StatDist = zeros(n_a,n_semiz,n_z,n_e,N_i,N_j);

%% Initialize distribution at age 1

for theta_c = 1:N_i
    % State variables, given PT and age j=1: (a,semiz,z,e)
    StatDist(:,:,:,:,theta_c,1) = joneDist*theta_dist(theta_c);
end

%% Forward induction
for theta_c = 1:N_i
    if simoptions.verbose==1; fprintf('Type %d out of %d \n',theta_c,N_i); end
    for j_c = 1:N_j-1
        if simoptions.verbose==1; fprintf('Age %d out of %d \n',j_c,N_j); end
        
        % Transition matrix for semiz is theta-dependent and age-dependent
        pi_semiz_ = pi_semiz_J(:,:,:,theta_c,j_c); % (semiz,semiz',f)
                      %1:3,a,semiz,z,e
        Policy_ = Policy(:,:,:,:,:,theta_c,j_c);
        [StatDist(:,:,:,:,theta_c,j_c+1)] = mu_onestep(StatDist(:,:,:,:,theta_c,j_c),Policy_,pi_semiz_,pi_z,pi_e);
    end %end loop over age
end % end loop over PT

%% Multiply with age distribution
for j_c = 1:N_j
     % (a,semiz,z,e,theta,j)
    StatDist(:,:,:,:,:,j_c) = StatDist(:,:,:,:,:,j_c)*mewj(j_c);
end

check_sum = sum(StatDist,"all");

if abs(check_sum-1)>1e-10
    warning('Distribution does not sum to one!')
end


end %end function fun_StatDist_cpu