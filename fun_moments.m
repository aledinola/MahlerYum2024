function [empl_rate_age] = fun_moments(mu,d1_grid,semiz_grid,pol_n,n_a,n_semiz,n_z,N_j,Names_i)

% (n,f,aprime,a,h,z,age,PT)
% States: a,h,z,age,PT
N_i = numel(Names_i);

n_grid = d1_grid;    % Grid for employment,    n=0,0.5,1.0
h_grid = semiz_grid; % Grid for health status: h=0,1,2

pol_n_val = n_grid(pol_n);

mass_alive_age = zeros(N_j,1);
for jj=1:N_j
    mass_alive_age(jj) = sum(mu(:,1:2,:,jj,:),"all");
end


empl_rate_age = zeros(N_j,1);

for jj=1:N_j
    empl_rate = 0;
    for ii=1:N_i
    for z_c=1:n_z
    for h_c=1:n_semiz
    for a_c=1:n_a
        n_val = pol_n_val(a_c,h_c,z_c,jj,ii);
        h_val = h_grid(h_c);
        alive = (h_val==0 || h_val==1);
        empl_rate = empl_rate+(n_val>0.0)*alive*mu(a_c,h_c,z_c,jj,ii);
    end
    end
    end
    end
    empl_rate_age(jj) = empl_rate;
end

empl_rate_age = empl_rate_age./mass_alive_age;

end %end function