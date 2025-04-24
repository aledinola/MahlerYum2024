function F = MahlerYum2024_ReturnFn(n,f,aprime,a,h,z,agej,Jr,lambda,wp,theta,...
    y_pen,r,kappa,nu_h0,nu_h1,iota_h0,iota_h1,gamma,psi,sigma,nu_e,b,educ)
%INPUT ARGUMENTS
% n:      Labor supply,             decision variable d1
% f:      health effort,            decision variable d2
% aprime: Next-period assets,       endogenous state choice
% a:      Current assets,           endogenous state
% h:      Health status,            semi-exogenous state
% z:      Labor productivity shock, exogenous state (Markov)
% agej:   Age
% lambda  Deterministic wage profile, depends on age and on PT educ=0,1
% wp:     Wage penalty, depends on PT educ=0,1
% nu_e:   Extra disutility of work if educ=0 (non-college), scalar
% educ:   Dummy 0,1 depends on PT educ=0,1

F=-Inf;

% Dead=zero utility
if h==2
    F=0;
    return
end

if h==0 
    % h=0, Unhealthy
    wage_penalty  = wp;      % Wage penalty
    kappa_shifter = kappa;   % Shifter for utility of consumption
    nu_shift      = nu_h0;   % Shifter for disutility of work
    iota_shift    = iota_h0; % Shifter for cost of effort
else
    % h=1, Healthy
    wage_penalty  = 0;       % Wage penalty
    kappa_shifter = 1;       % Shifter for utility of consumption
    nu_shift      = nu_h1;   % Shifter for disutility of work
    iota_shift    = iota_h1; % Shifter for cost of effort
end

if agej<Jr % If working age
    wage = exp(lambda*(1-wage_penalty)+theta+z);
    c    = (1+r)*a + wage*n-aprime;
else % Retirement
    c    = y_pen+(1+r)*a-aprime;
end

if c>0 
    work_disutil = nu_shift*exp(nu_e*(educ==0))*n^(1+1/gamma)/(1+1/gamma);
    effort_cost = iota_shift*f^(1+1/psi)/(1+1/psi);
    F= kappa_shifter*(c^(1-sigma)/(1-sigma)+b)-work_disutil-effort_cost;
end

end %end function