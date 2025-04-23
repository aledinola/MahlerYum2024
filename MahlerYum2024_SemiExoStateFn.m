function prob=MahlerYum2024_SemiExoStateFn(h,hprime,f,educ,eta,pi0,lambda1,delta,gamma1,gamma2)
%-------------------------------------------------------------------------%
% INPUT ARGUMENTS
% h:      Health status today, semi-exogenous state
% hprime: Health status tomorrow, semi-exogenous state
% f:      Health effort, decision variable d2
%-------------------------------------------------------------------------%
prob = -inf; % Just a placeholder (one that will cause errors if not overwritten)

% By how much exerting effort reduces the prob of bad health. If no effort,
% this is 1.
prob_bad_to_good = 1/(1+exp(-(pi0+lambda1*f+delta*h+gamma1*educ+gamma2*eta))  );
prob_good_to_good = prob_bad_to_good;

if h==0
    % Bad health today
    if hprime==0
        % bad health tomorrow
        prob = 1-prob_bad_to_good;
    elseif hprime==1
        % good health tomorrow
        prob = prob_bad_to_good;
    end
elseif h==1
    % Good health today
    if hprime==0
        % bad health tomorrow
        prob = 1-prob_good_to_good;
    elseif hprime==1
        prob = prob_good_to_good;
    end
end %end if zh

end %end function
