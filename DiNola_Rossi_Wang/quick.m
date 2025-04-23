% clear
% clc
% close all
% 
% pi0=1;
% pi1 = 0.5;
% fun = @(f) pi0*exp(-pi1*f);
% 
% f_grid = 0:0.01:1;
% figure
% plot(f_grid,fun(f_grid))

figure
plot(1:Params.J,mom.ave_age.frac_badhealth)
title('Fraction of sick people')

figure
plot(1:Params.J-1,mom.ave_age_health.effort(1,1:Params.J-1))
hold on
plot(1:Params.J-1,mom.ave_age_health.effort(2,1:Params.J-1))
legend('sick','healthy')
title('Effort by health')
