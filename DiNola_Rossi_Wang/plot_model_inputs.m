%% This script plots model inputs used later on in the computation

% Medical spending m(zh,e,j)
med_vec = zeros(G.n_z(2),G.n_e,G.N_j);
for jj=1:G.N_j
    med_vec(1,:,jj) = Params.medspend_j(jj)*G.e_grid;
end

figure
plot(Params.age,Params.kappa_j,'LineWidth',2)
xlabel('Age')
title('Age-dependent labor productivity')
grid on
print(fullfile(flag.results,'kappa_j'),'-dpng')

figure
plot(Params.age,Params.pen_j.low,'LineWidth',2)
hold on
plot(Params.age,Params.pen_j.high,'LineWidth',2)
legend('No college','College')
xlabel('Age')
ylabel('Monetary value')
title('Age-dependent pension benefits')
grid on
print(fullfile(flag.results,'pen_j'),'-dpng')

figure
for e_c = 1:G.n_e
    plot(Params.age,squeeze(med_vec(1,e_c,:)),'LineWidth',2)
    hold on
end
legend('e=1','e=2','e=3','e=4','e=5')
xlabel('Age')
ylabel('Monetary value')
title('Age-dependent medical expenses')
grid on

figure
plot(Params.age,G.zh_prob.low(1,:),'r','LineWidth',2)
hold on
plot(Params.age,G.zh_prob.high(1,:),'--r','LineWidth',2)
legend('Low \theta','high \theta','Location','northwest')
title('Share of Bad health households, for each \theta')
xlabel('Age')
ylabel('Share')
grid on

% Recall that zh_prob.PT has size (n_zh,N_age) and n_zh=2 (1=bad,2=good)
figure
plot(Params.age,G.zh_prob.low(1,:),'r','LineWidth',2)
hold on
plot(Params.age,G.zh_prob.low(2,:),'b','LineWidth',2)
hold on
plot(Params.age,G.zh_prob.high(1,:),'--r','LineWidth',2)
hold on
plot(Params.age,G.zh_prob.high(2,:),'--b','LineWidth',2)
legend('Low \theta - Bad health','low \theta - Good health',...
    'High \theta - Bad health','High \theta - Good health',...
    'Location', 'southoutside','NumColumns', 2)
title('Share of Good vs Bad health households, \theta=low')
xlabel('Age')
ylabel('Share')
grid on
print(fullfile(flag.results,'frac_health'),'-dpng')