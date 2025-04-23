function [prob] = Mod_SemiExoStateFn(zh,zh_prime,d2_effort,p_bb_j,p_gb_j,pi0,pi1)
% zh is *SEMI-EXOGENOUS* and *age-dependent*
% ------------------------------------------
% INPUTS
% zh:       Semi-exogenous health status today
% zh_prime: Semi-exogenous health status tomorrow
% d2_effort: Decision variable, effort
% OUTPUT
% prob:     Prob(zh_prime|zh), given effort d2
% ------------------------------------------

prob = -inf; % Just a placeholder (one that will cause errors if not overwritten)

% By how much exerting effort reduces the prob of bad health. If no effort,
% this is 1.
decrease = pi0*exp(-pi1*d2_effort);

if zh==0
    % Bad health today
    if zh_prime==0
        % bad health tomorrow
        prob = p_bb_j*decrease;
    elseif zh_prime==1
        % good health tomorrow
        prob = 1-p_bb_j*decrease;
    end
elseif zh==1
    % Good health today
    if zh_prime==0
        % bad health tomorrow
        prob = p_gb_j*decrease;
    elseif zh_prime==1
        prob = 1-p_gb_j*decrease;
    end
end %end if zh

end %end function