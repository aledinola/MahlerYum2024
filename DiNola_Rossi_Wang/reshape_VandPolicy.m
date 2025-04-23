function [V,mu,Policy] = reshape_VandPolicy(V_in,StatDist_in,Policy_in,n_a,n_semiz,n_z,n_e,N_i,N_j)
% Reshape stationary distribution and policy functions, to convert from
% toolkit notation (where PT are fields of struct) to better notation
%% State variables: (a,semiz,z,e,N_j) for each PT
% Assumptions: 
% only 2 permanent types N_i = (low and high), 
% 1 semiz shock
% 1 z shock, labor productivity

if iscell(N_i)
    N_i = numel(N_i);
end
if ~isstruct(Policy_in)
    error('Policy_in NOT a struct: reshape_VandPolicy assumes that you are using PTypes')
end
% 3 = choices are d1, d2 and a', in this order
if ~isequal(size(Policy_in.low),[3,n_a,n_semiz,n_z,n_e,N_j])
    error('Input variable Policy has wrong size!')
end

% Policy functions
Policy_low = squeeze(gather(Policy_in.low));
Policy_high = squeeze(gather(Policy_in.high));
%  I have to add extra dimension with N_i elements
Policy = zeros(3,n_a,n_semiz,n_z,n_e,N_i,N_j);
for j_c=1:N_j
    Policy(:,:,:,:,:,1,j_c) = Policy_low(:,:,:,:,:,j_c);
    Policy(:,:,:,:,:,2,j_c) = Policy_high(:,:,:,:,:,j_c);
end

% Value function
if ~isempty(V_in)
    if ~isstruct(V_in)
        error('V_in NOT a struct: reshape_VandPolicy assumes that you are using PTypes')
    end
    V_low = squeeze(gather(V_in.low));
    V_high = squeeze(gather(V_in.high));

    V = zeros(n_a,n_semiz,n_z,n_e,N_i,N_j);
    V(:,:,:,:,1,:) = V_low;
    V(:,:,:,:,2,:) = V_high;
else
    V = [];
end

% Distribution
if ~isempty(StatDist_in)
    if ~isstruct(StatDist_in)
        error('StatDist_in NOT a struct: reshape_VandPolicy assumes that you are using PTypes')
    end

    if ~isequal(size(StatDist_in.low),[n_a,n_semiz,n_z,n_e,N_j])
        error('Input variable StatDist has wrong size!')
    end

    StatDist_low = gather(StatDist_in.low);
    StatDist_high = gather(StatDist_in.high);

    mu = zeros(n_a,n_semiz,n_z,n_e,N_i,N_j);
    mu(:,:,:,:,1,:) = StatDist_low*StatDist_in.ptweights(1);
    mu(:,:,:,:,2,:) = StatDist_high*StatDist_in.ptweights(2);
else
    mu = [];
end

end %end function reshape_VandPolicy