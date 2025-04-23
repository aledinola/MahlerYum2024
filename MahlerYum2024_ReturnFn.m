function F = MahlerYum2024_ReturnFn(n,f,aprime,a,h,z,agej,Jr,lambda,wp,theta,...
    y_pen,r,kappa,nu_h0,nu_h1,iota_h0,iota_h1,gamma,psi,sigma,nu_e,b)
%INPUT ARGUMENTS
% n: Labor supply, decision variable d1
% f: health effort, decision variable d2
% aprime: Next-period assets, choice
% a: Current assets, endogenous state
% h: Health status, semi-exogenous state
% z: Labor productivity shock, exogenous state

F=-Inf;

if agej<Jr % If working age
    wage = exp(lambda*(1-wp*(h==0))+theta+z);
    c    = (1+r)*a + wage*n-aprime;
else % Retirement
    c    = y_pen+(1+r)*a-aprime;
end

% Shifter for utility of consumption
if h==0
    % Unhealthy
    kappa_shifter = kappa;
else
    kappa_shifter = 1;
end

% Shifter for disutility of work
if h==0
    nu_shift = nu_h0;
else
    nu_shift = nu_h1;
end

% Shifter for cost of effort
if h==0
    iota_shift = iota_h0;
else
    iota_shift = iota_h1;
end

if c>0 
    work_disutil = nu_shift*exp(nu_e)*n^(1+1/gamma)/(1+1/gamma);
    effort_cost = iota_shift*f^(1+1/psi)/(1+1/psi);
    F= kappa_shifter*(c^(1-sigma)/(1-sigma)+b)-work_disutil-effort_cost;
    % F=log(c/consumption_equiv_units)+psi*log(leisure)+utility_of_children;
end

end %end function