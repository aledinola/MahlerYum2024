function [] = plot_model_results(Params,G,mom,flag)

%  Unpack
N_j = G.N_j;
ave            = mom.ave;
ave_age        = mom.ave_age;
%ave_health     = mom.ave_health;
ave_age_health = mom.ave_age_health;


var_names = {'frac_badhealth','cons','labincome','hours','assets','medical'};

% Moments by age 
for ii=1:numel(var_names)
    figure
    name_ii = var_names{ii};
    plot(1:N_j,ave_age.(name_ii))
    yline(ave.(name_ii))
    legend('By age','Unconditional mean')
    title(name_ii)
    xlabel('Age, j')
    ylabel('Mean')
    fig_name = [name_ii,'_age'];
    print(fullfile(flag.results,fig_name),'-dpng')
end

% Moments by age and health
for ii=2:numel(var_names)
    figure
    name_ii = var_names{ii};
    plot(1:N_j,ave_age_health.(name_ii)(1,:))
    hold on
    plot(1:N_j,ave_age_health.(name_ii)(2,:))
    title(name_ii)
    legend('Sick','Healthy')
    xlabel('Age, j')
    ylabel('Mean')
    fig_name = [name_ii,'_age_health'];
    print(fullfile(flag.results,fig_name),'-dpng')
end

end %end function
