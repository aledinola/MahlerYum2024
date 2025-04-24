function prob=MahlerYum2024_SemiExoStateFn(h,hprime,f,sj_h0,sj_h1,educ,eta,pi0,lambda1,delta,gamma1,gamma2)
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS
% h:      Health status today, semi-exogenous state
%         0=unhealthy, 1=healthy, 2=dead
% hprime: Health status tomorrow, semi-exogenous state
% f:      Health effort, decision variable d2
% sj_h0:  Survival probability, given h=0 (unhealthy), depends on educ and age
% sj_h1:  Survival probability, given h=1 (healthy), depends on educ and age
% educ:   Permanent education type, as dummy 0,1
% eta:    Permanent healthy type, as dummy 0,1
% pi0:
% lambda1,delta,gamma1,gamma2: Scalar parameters
%-------------------------------------------------------------------------%
prob = -inf; % Just a placeholder (one that will cause errors if not overwritten)

% Probability of being healthy in the next period, Pr(hprime=1|h,f,...)
p_hprime_1 = 1/(1+exp(-(pi0+lambda1*f+delta*h+gamma1*educ+gamma2*eta))  );

if h==0
    % Unhealthy current period
    if hprime==0
        % Unhealthy next period
        prob = sj_h0*(1-p_hprime_1);
    elseif hprime==1
        % Healthy next period
        prob = sj_h0*p_hprime_1;
    elseif hprime==2
        % Dead next period
        prob = 1-sj_h0;
    end
elseif h==1
    % Healthy current period
    if hprime==0
        % Unhealthy next period
        prob = sj_h1*(1-p_hprime_1);
    elseif hprime==1
        % Healthy next period
        prob = sj_h1*p_hprime_1;
    elseif hprime==2
        % Dead next period
        prob = 1-sj_h1;
    end
elseif h==2
    % Dead in current period: absorbing state
    if hprime==0
        prob = 0;
    elseif hprime==1
        prob = 0;
    elseif hprime==2
        prob = 1;
    end
end %end if zh

end %end function
