function labincome = labincome_fn(n,f,aprime,a,h,z,agej,Jr,lambda_h0,lambda_h1,theta)
%INPUT ARGUMENTS
% n:      Labor supply,             decision variable d1
% f:      health effort,            decision variable d2
% aprime: Next-period assets,       endogenous state choice
% a:      Current assets,           endogenous state
% h:      Health status,            semi-exogenous state
% z:      Labor productivity shock, exogenous state (Markov)
% --- Parameters
% VARIABLE   DESCRIPTION        AGE-DEPENDENT      TYPE-DEPENDENT
% agej:      Age                     Yes                No
% Jr:        Retirement age          No                 No
% lambda_h0  Deterministic wage,     Yes                Yes (e=0,1)
% lambda_h1  Deterministic wage,     Yes                Yes (e=0,1)
% nu_e:      Extra disut. of work    No                 No
% educ:      Dummy 0,1               No                 Yes (e=0,1)


% Dead=zero utility
if h==2
    labincome=0;
    return
end

if h==0 
    % h=0, Unhealthy
    lambda = lambda_h0;   % Deterministic component of earnings
else
    % h=1, Healthy
    lambda = lambda_h1; % Deterministic component of earnings
end

if agej<Jr % If working age
    wage = exp(lambda+theta+z);
    labincome = wage*n;
else % Retirement
    labincome = 0;
end

end %end function