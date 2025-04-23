function val = fun_values_on_grid(Params,Policy_in,FnsToEval,d_grid,a_grid,semiz_grid,z_grid,e_grid,n_d,n_a,n_semiz,n_z,n_e,N_i,N_j)

if iscell(N_i)
    n_theta = numel(N_i);
else
    n_theta = N_i;
end
theta_grid = Params.theta_i;

d1_grid = d_grid(1:n_d(1)); %Grid for labor supply, n
d2_grid = d_grid(n_d(1)+1:sum(n_d)); %Grid for health effort, f

% Policy functions: reshape them but keep them on the GPU
Policy_low = squeeze(Policy_in.low);
Policy_high = squeeze(Policy_in.high);
%  I have to add extra dimension with N_i elements
Policy = zeros(3,n_a,n_semiz,n_z,n_e,n_theta,N_j,'gpuArray');
for j_c=1:N_j
    Policy(:,:,:,:,:,1,j_c) = Policy_low(:,:,:,:,:,j_c);
    Policy(:,:,:,:,:,2,j_c) = Policy_high(:,:,:,:,:,j_c);
end

% Policy for d1, d2 and assets, dim: (a,semiz,z,e,theta_i,age)
pol_d1     = reshape(Policy(1,:,:,:,:,:,:),n_a,n_semiz,n_z,n_e,n_theta,N_j);
pol_d2     = reshape(Policy(2,:,:,:,:,:,:),n_a,n_semiz,n_z,n_e,n_theta,N_j);
pol_aprime = reshape(Policy(3,:,:,:,:,:,:),n_a,n_semiz,n_z,n_e,n_theta,N_j);

% Now apply grid: indices-->values
pol_d1 = gpuArray(d1_grid(pol_d1));
pol_d2 = gpuArray(d2_grid(pol_d2));
pol_aprime = gpuArray(a_grid(pol_aprime));

a      = gpuArray(a_grid);                   %(a,1,    1,1)
semiz  = gpuArray(shiftdim(semiz_grid,-1));  %(1,semiz,1,1)
z      = gpuArray(shiftdim(z_grid,-2));      %(1,1,    z,1)
e      = gpuArray(shiftdim(e_grid,-3));      %(1,1,    1,e)

% Loop over ages and over type
names = fieldnames(FnsToEval);

for ii=1:numel(names)
    name_var = names{ii};
    val.(name_var) = zeros(n_a,n_semiz,n_z,n_e,n_theta,N_j);
    for j_c=1:N_j
        for theta_c=1:n_theta
            theta_i = theta_grid(theta_c);

            % Extract all parameters by PT and by age
            % kappa_j = Params.kappa_j(j_c);
            % r       = Params.r;
            % w       = Params.w;
            P = fun_params_cpu(Params,theta_c,j_c);

            d1     = pol_d1(:,:,:,:,theta_c,j_c);     %(a,semiz,z,e)
            d2     = pol_d2(:,:,:,:,theta_c,j_c);     %(a,semiz,z,e)
            aprime = pol_aprime(:,:,:,:,theta_c,j_c); %(a,semiz,z,e)

            val.(name_var)(:,:,:,:,theta_c,j_c) = arrayfun(FnsToEval.(name_var),d1,d2,aprime,...
                a,semiz,z,e,theta_i,P.kappa_j,P.pen_j,P.medspend_j,P.varrho,P.w,P.r,P.c_floor,...
                P.deduc,P.tau0,P.tau1,P.sub,P.tau_ss,P.cap_ss);
        end % end PT
    end % end age j
end %end names

end %end function