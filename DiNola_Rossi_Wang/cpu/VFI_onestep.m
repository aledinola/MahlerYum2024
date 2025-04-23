function [V,Policy] = VFI_onestep(V_next,semiz_grid,pi_semiz,a_grid,z_grid,pi_z,e_grid,pi_e,d_grid,theta,P,n_d)
% PURPOSE
%  Onestep iteration of Bellman equation, model with semi-exogenous shock
% INPUTS
%  V_next(a',semiz',z',e')
%  semiz_grid is (n_semiz,1)
%  pi_semiz is (semiz,semiz',f)

% Initialize output
%V = zeros(n_a,n_semiz,n_z,n_e);
%Policy = zeros(2,n_a,n_semiz,n_z,n_e);

if ~isstruct(P)
    error('VFI_onestep: Input P must be a structure')
end

% Discount factor and age-dependent survival prob
beta = P.beta;
sj   = P.sj;

N_d  = size(d_grid,1);
n_a  = size(a_grid,1);
n_semiz = size(semiz_grid,1);
n_z     = size(z_grid,1);
n_e  = size(e_grid,1);

%N_d = n_d(1)*n_d(2);

aprime = a_grid;
a      = a_grid';
V_d      = zeros(n_a,n_semiz,n_z,n_e,N_d);
Pol_aprime_d = zeros(n_a,n_semiz,n_z,n_e,N_d); % assets, given d

for d_c=1:N_d
    d1 = d_grid(d_c,1); % labor supply, n
    d2 = d_grid(d_c,2); % health effort, f
    [~,d2_c] = ind2sub([n_d(1),n_d(2)],d_c);
    pi_semiz_d = pi_semiz(:,:,d2_c);
    for e_c=1:n_e
    e = e_grid(e_c);
    for z_c=1:n_z
    z = z_grid(z_c);
    for semiz_c=1:n_semiz
        semiz = semiz_grid(semiz_c);
        % EVz is (a',1)
        EVz = compute_EV(V_next,pi_semiz_d(semiz_c,:),pi_z(z_c,:),pi_e,n_a,n_semiz,n_z,n_e);
        % RetMat is (a',a)
        RetMat = Mod_ReturnFn_cpu(d1,d2,aprime,a,semiz,z,e,theta,P.kappa_j,...
            P.pen_j,P.medspend_j,P.varrho,P.w,P.r,P.c_floor,P.deduc,P.tau0,...
            P.tau1,P.sub,P.tau_ss,P.cap_ss,P.sigma,P.delta,P.nu,P.gamma,P.iota_j,P.psi);
        %entire_rhs is (a',a) since RetMat is (a',a)
        entire_rhs = RetMat+beta*sj*EVz;
        [max_val,max_ind] = max(entire_rhs,[],1); %(1,a)
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

end %end function "VFI_onestep"

function EVz = compute_EV(V_next,prob1,prob2,prob3,n_a,nz1,nz2,nz3)
% V_next is (a',z1',z2',z3')
% Assumption: prob1,prob2,prob3 are vectors (row or column does not matter)

EVz = zeros(n_a,1);

for z3_c = 1:nz3
    for z2_c = 1:nz2
        for z1_c = 1:nz1
            EVz = EVz + V_next(:,z1_c,z2_c,z3_c)*prob1(z1_c)*prob2(z2_c)*prob3(z3_c);
        end
    end
end

end %end subfunction "compute_EV"