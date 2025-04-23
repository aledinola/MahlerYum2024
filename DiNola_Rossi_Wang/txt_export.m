function [] = txt_export(Params,Gr,mom,flag)
% In this function we export model targets in txt files. Also display on
% screen
% State variables:
% (a,z1,z2,e,PT,j)

if ~isstruct(Params)
    error('Input Params must be a structure')
end
if ~isstruct(Gr)
    error('Input Gr must be a structure')
end
if ~isstruct(mom)
    error('Input mom must be a structure')
end
if ~isstruct(flag)
    error('Input flag must be a structure')
end

% Unpack variables from struct "mom"
ave = mom.ave;
ave_young = mom.ave_young;
ave_health = mom.ave_health;
ave_health_young = mom.ave_health_young;
gini = mom.gini;
gini_health = mom.gini_health;
corr = mom.corr;

fid=fopen(fullfile(flag.results,'targets_model_manual.txt'),'wt');  % overwrite

fprintf(fid,' \n');
fprintf(fid,'GRIDS SIZES \n');
fprintf(fid,'n_d1,  hours of work:   %d \n',Gr.n_d(1));
fprintf(fid,'n_d2,  health effort:   %d \n',Gr.n_d(2));
fprintf(fid,'n_a,  assets:          %d \n',Gr.n_a); 
fprintf(fid,'n_semiz, health shock:    %d \n',Gr.n_semiz);
fprintf(fid,'n_z, prod shock:      %d \n',Gr.n_z);
fprintf(fid,'n_e,  medical shock:   %d \n',Gr.n_e);
fprintf(fid,'J,    periods in life: %d \n',Gr.N_j);

fprintf(fid,'PARAMETERS \n');
fprintf(fid,'Discount factor, beta:        %f \n',Params.beta); 
fprintf(fid,'CRRA consumption, sigma:      %f \n',Params.sigma); 
fprintf(fid,'Elast labor supply, gamma: %f \n',Params.gamma);
fprintf(fid,'Weight on labor disutil, nu:  %f \n',Params.nu);
fprintf(fid,'Elast of health effort, psi: %f \n',Params.psi);
fprintf(fid,'Weight on health effort disutil, iota, j=1:  %f \n',Params.iota_j(1));
fprintf(fid,'Weight on health effort disutil, iota, j=J:  %f \n',Params.iota_j(Params.J));
fprintf(fid,'Health transition, pi0: %f \n',Params.pi0);
fprintf(fid,'Health transition, pi1: %f \n',Params.pi1);

fprintf(fid,'Prob sick to sick, p_bb, j=1:    %f \n',Params.p_bb_j(1));
fprintf(fid,'Prob sick to sick, p_bb, j=J:     %f \n',Params.p_bb_j(Params.J));
fprintf(fid,'Prob healthy to sick, p_gb, j=1:  %f \n',Params.p_gb_j(1,1));
fprintf(fid,'Prob healthy to sick, p_gb, j=J:  %f \n',Params.p_gb_j(1,Params.J));

fprintf(fid,'Decrease in product if sick, low PT: %f \n',Params.varrho.low);
fprintf(fid,'Decrease in product if sick, high PT: %f \n',Params.varrho.high);
fprintf(fid,'Preference shifter delta: %f \n',Params.delta);

fprintf(fid,'Tax level tau0:         %f \n',Params.tau0); 
fprintf(fid,'Tax progressivity tau1: %f \n',Params.tau1); 
fprintf(fid,'Subsidy medical expens: %f \n',Params.sub);
fprintf(fid,'Soc secur tax rate:     %f \n',Params.tau_ss); 
fprintf(fid,'Soc secur income cap:   %f \n',Params.cap_ss); 
fprintf(fid,'Consump floor, c_floor: %f \n',Params.c_floor); 

fprintf(fid,'RESULTS \n');
fprintf(fid,'--- Government ---- \n');
fprintf(fid,'atr of mean income:    %f \n',mom.tax.atr_mean); 
fprintf(fid,'Gov net tax revenue:   %f \n',mom.tax.govnet); 
fprintf(fid,'--- Unconditional Means ---- \n');
fprintf(fid,'Aggregate welfare:     %f \n',ave.welfare);
fprintf(fid,'Average assets:        %f \n',ave.assets); 
fprintf(fid,'Average hours (non-retirees): %f \n',ave_young.hours);
fprintf(fid,'Average effort: %f \n',ave.effort);
fprintf(fid,'Average labincome:     %f \n',ave.labincome);
fprintf(fid,'Average income:        %f \n',ave.income);
fprintf(fid,'Average consumption:   %f \n',ave.cons);
fprintf(fid,'Total tax revenues:    %f \n',ave.taxrev);
fprintf(fid,'Average share sick:    %f \n',ave.frac_badhealth);
fprintf(fid,'Total medical expens:  %f \n',ave.medical);
fprintf(fid,'Welfare transfers b:   %f \n',ave.b_tr);

fprintf(fid,'--- Means conditional on health ---- \n');
fprintf(fid,'Average assets, sick:       %f \n',ave_health.assets(1)); 
fprintf(fid,'Average assets, healthy:    %f \n',ave_health.assets(2)); 
fprintf(fid,'Average hours, sick (non-retirees):        %f \n',ave_health_young.hours(1)); 
fprintf(fid,'Average hours, healthy (non-retirees):     %f \n',ave_health_young.hours(2)); 
fprintf(fid,'Average health effort, sick:       %f \n',ave_health.effort(1)); 
fprintf(fid,'Average health effort, healthy:    %f \n',ave_health.effort(2)); 
fprintf(fid,'Average labincome, sick:    %f \n',ave_health.labincome(1)); 
fprintf(fid,'Average labincome, healthy: %f \n',ave_health.labincome(2));
fprintf(fid,'Average income, sick:       %f \n',ave_health.income(1)); 
fprintf(fid,'Average income, healthy:    %f \n',ave_health.income(2));
fprintf(fid,'Average cons, sick:         %f \n',ave_health.cons(1)); 
fprintf(fid,'Average cons, healthy:      %f \n',ave_health.cons(2));
fprintf(fid,'Average medical expens, sick:    %f \n',ave_health.medical(1)); 
fprintf(fid,'Average medical expens, healthy: %f \n',ave_health.medical(2));
fprintf(fid,'Welfare transfers b, sick:       %f \n',ave_health.b_tr(1));
fprintf(fid,'Welfare transfers b, healthy:    %f \n',ave_health.b_tr(2));

fprintf(fid,'--- Correlations, unconditional ---- \n');
fprintf(fid,'labincome and health:        %f \n',corr.labincome_health); 

fprintf(fid,'--- Unconditional Gini ---- \n');
fprintf(fid,'Gini assets:        %f \n',gini.assets); 
%fprintf(fid,'Gini hours:         %f \n',gini.hours);
fprintf(fid,'Gini effort:         %f \n',gini.effort);
fprintf(fid,'Gini labincome:     %f \n',gini.labincome);
fprintf(fid,'Gini income:        %f \n',gini.income);
fprintf(fid,'Gini consumption:   %f \n',gini.cons);
fprintf(fid,'Gini tax revenues:    %f \n',gini.taxrev);
fprintf(fid,'Gini medical expens:  %f \n',gini.medical);

fprintf(fid,'--- Gini conditional on health ---- \n');
fprintf(fid,'Gini assets, sick:       %f \n',gini_health.assets(1)); 
fprintf(fid,'Gini assets, healthy:    %f \n',gini_health.assets(2)); 
fprintf(fid,'Gini effort, sick:        %f \n',gini_health.effort(1)); 
fprintf(fid,'Gini effort, healthy:     %f \n',gini_health.effort(2)); 
fprintf(fid,'Gini labincome, sick:    %f \n',gini_health.labincome(1)); 
fprintf(fid,'Gini labincome, healthy: %f \n',gini_health.labincome(2));
fprintf(fid,'Gini income, sick:       %f \n',gini_health.income(1)); 
fprintf(fid,'Gini income, healthy:    %f \n',gini_health.income(2));
fprintf(fid,'Gini cons, sick:         %f \n',gini_health.cons(1)); 
fprintf(fid,'Gini cons, healthy:      %f \n',gini_health.cons(2));
fprintf(fid,'Gini medical expens, sick:    %f \n',gini_health.medical(1)); 
fprintf(fid,'Gini medical expens, healthy: %f \n',gini_health.medical(2));
fclose(fid);

end %end function