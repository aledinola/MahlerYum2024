function [V,Policy] = fun_vfi_cpu(Params,Gr,vfoptions)

% --- Grid dimensions
n_d     = Gr.n_d;
n_a     = Gr.n_a;
n_z     = Gr.n_z;          % No. grid points for z, exogenous shock
n_e     = Gr.n_e;
n_theta = Gr.n_theta;
N_j     = Gr.N_j;
n_semiz = Gr.n_semiz;  % No. grid points for semiz, semi-exogenous shock
% --- Grids
d_grid_in  = Gr.d_grid;
a_grid     = Gr.a_grid; % Grid for assets
theta_grid = Params.theta_i;
semiz_grid = vfoptions.semiz_grid; % (n_semiz,1)
z_grid     = Gr.z_grid;     % (n_z1,1)
e_grid     = Gr.e_grid;     % (n_e,1)
pi_z       = gather(Gr.pi_z);       % (n_z1,n_z1)
pi_e       = gather(Gr.pi_e);       % (n_e,1)
% --- The transition matrix depends on effort, type theta and age j
pi_semiz_J = gather(Gr.pi_semiz_J); % (semiz,semiz'|f,theta,j)

% Combine two decision variables (n,f) into d=[d1,d2] where only d2 affects
% pi_semiz
n_grid = d_grid_in(1:n_d(1));
f_grid = d_grid_in(n_d(1)+1:sum(n_d));
[n_grid1,f_grid1]=ndgrid(n_grid,f_grid);
n_grid1=n_grid1(:);
f_grid1=f_grid1(:);

d_grid = [n_grid1,f_grid1];
if ~isequal(size(d_grid),[n_d(1)*n_d(2),2])
    error('Homemade d_grid has wrong dimensions')
end

V      = zeros(n_a,n_semiz,n_z,n_e,n_theta,N_j);
% The first dim of Policy contains policies for d1,d2,a'
Policy = zeros(3,n_a,n_semiz,n_z,n_e,n_theta,N_j);

%% Last period of life

if vfoptions.verbose==1; fprintf('Age %d out of %d \n',N_j,N_j); end

for theta_c = 1:n_theta
    theta_val = theta_grid(theta_c);
    P = fun_params_cpu(Params,theta_c,N_j);
    [V_last,Policy_last] = solve_last_period(a_grid,d_grid,semiz_grid,z_grid,e_grid,n_d,P,theta_val);
    V(:,:,:,:,theta_c,N_j) = V_last;
    Policy(:,:,:,:,:,theta_c,N_j) = Policy_last;
end %end theta

%% Backward iteration

for theta_c = 1:n_theta
    if vfoptions.verbose==1; fprintf('Type %d out of %d \n',theta_c,n_theta); end
    theta = theta_grid(theta_c);
    
    for j_c = N_j-1:-1:1

    if vfoptions.verbose==1; fprintf('Age %d out of %d \n',j_c,N_j); end

    % Extract parameters that depend on age j and type theta_c in P
    P = fun_params_cpu(Params,theta_c,j_c);

    % V_next(a',semiz',z',e')
    V_next = V(:,:,:,:,theta_c,j_c+1);

    % Transition matrix for semiz is theta-dependent and age-dependent
    pi_semiz_ = pi_semiz_J(:,:,:,theta_c,j_c); % (semiz,semiz',f)

    % Call one-step function
    [V(:,:,:,:,theta_c,j_c),Policy(:,:,:,:,:,theta_c,j_c)] = VFI_onestep(V_next,...
        semiz_grid,pi_semiz_,a_grid,z_grid,pi_z,e_grid,pi_e,d_grid,theta,P,n_d);

    end
end

% Policy(1,etc.) is for d
% Policy(2,etc.) is for a'

end %end function

function [V,Policy] = solve_last_period(a_grid,d_grid,semiz_grid,z_grid,e_grid,n_d,P,theta)

N_d     = size(d_grid,1);
n_a     = size(a_grid,1);
n_semiz = size(semiz_grid,1);
n_z     = size(z_grid,1);
n_e     = size(e_grid,1);

aprime = a_grid;
a      = a_grid';
V_d      = zeros(n_a,n_semiz,n_z,n_e,N_d);
Pol_aprime_d = zeros(n_a,n_semiz,n_z,n_e,N_d); % assets, given d

for d_c=1:N_d
    d1 = d_grid(d_c,1); % labor supply, n
    d2 = d_grid(d_c,2); % health effort, f
    for e_c=1:n_e
    e = e_grid(e_c);
    for z_c=1:n_z
    z = z_grid(z_c);
    for semiz_c=1:n_semiz
        semiz = semiz_grid(semiz_c);
        RetMat = Mod_ReturnFn_cpu(d1,d2,aprime,a,semiz,z,e,theta,P.kappa_j,...
            P.pen_j,P.medspend_j,P.varrho,P.w,P.r,P.c_floor,P.deduc,P.tau0,...
            P.tau1,P.sub,P.tau_ss,P.cap_ss,P.sigma,P.delta,P.nu,P.gamma,P.iota_j,P.psi);
        [max_val,max_ind] = max(RetMat,[],1); %(1,a)
        V_d(:,semiz_c,z_c,e_c,d_c) = max_val;          % best V, given d
        Pol_aprime_d(:,semiz_c,z_c,e_c,d_c) = max_ind; % best a', given d
    end %end semiz
    end %end z
    end %end e
end %end d

% Now maximize with respect to d
[V,d_max] = max(V_d,[],5); % this gives (a,semiz,z,e)
% 1st slot for d1, 2nd slot for d2, 3rd slot for a'
Policy = zeros(3,n_a,n_semiz,n_z,n_e);
% Size of d_max is (a,semiz,z,e)
[d1_max,d2_max] = ind2sub([n_d(1),n_d(2)],d_max);
Policy(1,:,:,:,:) = d1_max;
Policy(2,:,:,:,:) = d2_max;
for e_c=1:n_e
for z_c=1:n_z
for semiz_c=1:n_semiz
for a_c=1:n_a
    Policy(3,a_c,semiz_c,z_c,e_c) = Pol_aprime_d(a_c,semiz_c,z_c,e_c,d_max(a_c,semiz_c,z_c,e_c));
end
end
end
end

end %end subfunction