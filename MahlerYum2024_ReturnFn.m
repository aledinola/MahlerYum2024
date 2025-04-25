function F = MahlerYum2024_ReturnFn(n,f,aprime,a,h,z,agej,Jr,lambda_h0,lambda_h1,...
    theta,y_pen,r,kappa_tilde,nu_h0,nu_h1,iota_h0,iota_h1,gamma,psi,sigma,...
    nu_e,b,educ,tau_s,tau_p,ybar)
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

F=-Inf;

% Dead=zero utility
if h==2
    F=0;
    return
end

if h==0 
    % h=0, Unhealthy
    lambda        = lambda_h0;   % Deterministic component of earnings
    kappa_shifter = kappa_tilde; % Shifter for utility of consumption
    nu_shift      = nu_h0;       % Shifter for disutility of work
    iota_shift    = iota_h0;     % Shifter for cost of effort
else
    % h=1, Healthy
    lambda        = lambda_h1; % Deterministic component of earnings
    kappa_shifter = 1;         % Shifter for utility of consumption
    nu_shift      = nu_h1;     % Shifter for disutility of work
    iota_shift    = iota_h1;   % Shifter for cost of effort
end

if agej<Jr % If working age
    wage = exp(lambda+theta+z);
    earnings = wage*n;
    %taxab_income = earnings+r*a;
    %tax=taxab_income-(1-tau_s)*taxab_income^(1-tau_p) *(ybar^tau_p);
    c    = (1+r)*a + earnings - aprime;
else % Retirement
    c    = y_pen + (1+r)*a - aprime;
end

if c>0 
    util_cons    = kappa_shifter*(c^(1-sigma)/(1-sigma)+b);
    work_disutil = nu_shift*exp(nu_e*(educ==0))*n^(1+1/gamma)/(1+1/gamma);
    effort_cost  = iota_shift*f^(1+1/psi)/(1+1/psi);
    F= util_cons-work_disutil-effort_cost;
end

end %end function