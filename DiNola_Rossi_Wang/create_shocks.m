function [p_bb_j,p_gb_j,zh_grid,dist_health,z_grid,pi_z,e_grid,pi_e] = create_shocks(n_semiz,n_z,n_e,N_j,Par)
% PURPOSE
%  This function creates grids and transition matrices for income (n1 points) 
%  and health (n2 points) shocks. Note that health shock depends on age and 
%  on permanent type theta.
% OUTPUTS
%  p_bb_j,p_gb_j: Age-dependent parameters for health transition 
%  zh_grid:       Grid for health, semi-exog shock (trans probs are defined
%                 in function Mod_SemiExoStateFn)
%  dist_health:   Initial distribution of health at age j=1 
%  z_grid,pi_z:   Grid and transition matrix for exogen product shock zn
%  e_grid,pi_e:   Pure iid shock e for medical costs


zh_grid = linspace(0,1,n_semiz)';

%% Income shock
[z_grid,pi_z]=discretizeAR1_Rouwenhorst(0.0,Par.rho,Par.sigma_eps,n_z);
z_grid = exp(z_grid);

%% Health shock, age-dependent and Ptype-dependent {low,high}
% 0=bad health, 1=good health
dist_health = [0.05,0.95]; % 5% in bad health at age j=1
% NOTE: If the grid zh does not depend on age (only transition matrix
% does), should I still replicate along 2nd dimension?
% zh_grid_J.low  = repmat([0,1]',1,N_j); % [n2,N_j]
% zh_grid_J.high = repmat([0,1]',1,N_j); % [n2,N_j]
% % bb: Prob of having bad health tomorrow given that health today is bad
p_bb = linspace(0.6,0.9,N_j);
% % gb: Prob of having bad health tomorrow given that health today is good
% % 0.1-0.4 for theta(1) and 0.0-0.2 for theta(2)
p_gb_low  = linspace(0.1,0.4,N_j);
p_gb_high = linspace(0.1,0.4,N_j);
% pi_zh_J.low  = zeros(n2,n2,N_j);
% pi_zh_J.high = zeros(n2,n2,N_j);
% for jj=1:N_j
%     % -- Low type (no college educ)
%     prob_j = [p_bb(jj),1-p_bb(jj);
%               p_gb_low(jj),1-p_gb_low(jj)];
%     pi_zh_J.low(:,:,jj) = prob_j;
%     % -- High type (college educ)
%     prob_j = [p_bb(jj),1-p_bb(jj);
%              p_gb_high(jj),1-p_gb_high(jj)];
%     pi_zh_J.high(:,:,jj) = prob_j;
% end
% disp('Check pi_zh_J.low...')
% check_markov_age(pi_zh_J.low,n_z(2),N_j)
% disp('Check pi_zh_J.high...')
% check_markov_age(pi_zh_J.high,n_z(2),N_j)
% 
% % Recall two permanent types: low,high
% z_grid_J.low=[zn_grid*ones(1,N_j); zh_grid_J.low];   % array with size [n1+n2,N_j]
% z_grid_J.high=[zn_grid*ones(1,N_j); zh_grid_J.high]; % array with size [n1+n2,N_j]
% pi_z_J.low=zeros(n1*n2,n1*n2,N_j);  % array with size [n1*n2,n1*n2,N_j]
% pi_z_J.high=zeros(n1*n2,n1*n2,N_j); % array with size [n1*n2,n1*n2,N_j]
% for jj=1:N_j
%    pi_z_J.low(:,:,jj)=kron(pi_zh_J.low(:,:,jj), pi_zn);  % note: kron in reverse order
%    pi_z_J.high(:,:,jj)=kron(pi_zh_J.high(:,:,jj), pi_zn);  % note: kron in reverse order
% end
% disp('Check pi_z_J.low...')
% check_markov_age(pi_z_J.low,prod(n_z),N_j)
% disp('Check pi_z_J.high...')
% check_markov_age(pi_z_J.high,prod(n_z),N_j)
% 
% if ~isstruct(z_grid_J)
%     error('Output z_grid_J must be a structure')
% end
% if ~isstruct(pi_z_J)
%     error('Output pi_z_J must be a structure')
% end

%% Iid shock e
%Tauchen_q = 3.0;
%[e_grid,pi_e]=discretizeAR1_Tauchen(0.0,0.0,Par.sig_e,n_e,Tauchen_q);
%e_grid = exp(e_grid);
%pi_e = pi_e(1,:)';
[e_grid, pi_e] = normal_discrete_1(n_e, 0.0, Par.sig_e);
e_grid = exp(e_grid);
pi_e = pi_e/sum(pi_e);
%e_grid = [1,1]';
%pi_e = [1,0]';

% Pass also age-dependent
p_bb_j = p_bb; %dim: (1,N_j)
p_gb_j = [p_gb_low;p_gb_high]; %dim: (N_i,N_j)

end %end function create_shocks