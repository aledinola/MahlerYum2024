function [mu,pol_n,pol_f,pol_aprime] = create_my_statdist(StatDist,Policy,Names_i,n_a,n_semiz,n_z,N_j)

N_i = numel(Names_i);

ptweights = StatDist.ptweights;

% Create stationary distribution
mu = zeros(n_a,n_semiz,n_z,N_j,N_i);
for ii=1:N_i
    mu(:,:,:,:,ii) = StatDist.(Names_i{ii})*ptweights(ii);
end

% Create policy functions for n,f and aprime (indexes)
pol_n = zeros(n_a,n_semiz,n_z,N_j,N_i);
pol_f = zeros(n_a,n_semiz,n_z,N_j,N_i);
pol_aprime = zeros(n_a,n_semiz,n_z,N_j,N_i);

for ii=1:N_i
    pol_n(:,:,:,:,ii) = squeeze(Policy.(Names_i{ii})(1,:,:,:,:,:));
    pol_f(:,:,:,:,ii) = squeeze(Policy.(Names_i{ii})(2,:,:,:,:,:));
    pol_aprime(:,:,:,:,ii) = squeeze(Policy.(Names_i{ii})(3,:,:,:,:,:));
end

check_sum = sum(mu,"all");
if abs(check_sum-1)>1e-8
    error('Distribution mu does not sum to one!')
end

end %end function