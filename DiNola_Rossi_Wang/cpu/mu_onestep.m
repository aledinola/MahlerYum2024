function StatDist_new = mu_onestep(StatDist,Policy_,pi_semiz_,pi_z,pi_e)
% StatDist is dim (a,semiz,z,e)
% pi_semiz_ is (semiz,semiz',f)
% pi_z is (z,z')
% pi_e  is (e',1)

[n_a,n_semiz,n_z,n_e] = size(StatDist);

% Extract policies for d2 and a' from Policy_
Pol_d2     = reshape(Policy_(2,:,:,:,:),[n_a,n_semiz,n_z,n_e]); %index
Pol_aprime = reshape(Policy_(3,:,:,:,:),[n_a,n_semiz,n_z,n_e]); %index

% Initialize output
StatDist_new = zeros(n_a,n_semiz,n_z,n_e);

for e_c = 1:n_e
for z_c = 1:n_z
for semiz_c = 1:n_semiz
for a_c = 1:n_a
    ap_opt = Pol_aprime(a_c,semiz_c,z_c,e_c);
    d2_opt = Pol_d2(a_c,semiz_c,z_c,e_c);
    pi_semiz = pi_semiz_(:,:,d2_opt);
    for ep_c = 1:n_e
    for zp_c = 1:n_z
    for semizp_c = 1:n_semiz
        StatDist_new(ap_opt,semizp_c,zp_c,ep_c) = StatDist_new(ap_opt,semizp_c,zp_c,ep_c)...
            +pi_semiz(semiz_c,semizp_c)*pi_z(z_c,zp_c)*pi_e(ep_c)*StatDist(a_c,semiz_c,z_c,e_c);
    end
    end
    end
end %end a
end %end semiz
end %end z
end %end e
% (semiz,z,a)-->(a,semiz,z)
%StatDist_new = permute(StatDist_new,[3,1,2]);

check = sum(StatDist_new,"all");

end %end function