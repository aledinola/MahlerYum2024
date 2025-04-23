function cons=Mod_cons(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,...
    varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss)
% INPUTS
%  n:       Labor supply (full-time equivalents),d1
%  f:       Health effort, d2
%  aprime:  Next-period assets
%  a:       Current assets
%  zh:      Health shock (0=sick, 1=healthy)
%  zn:      Labor productivity shock
%  e:       Transitory iid medical shock
%  theta_i: Permanent type (education)
% OUTPUT
%  cons:    Consumption

% If zh=0 (bad health), there is an income penalty. Note that varrho<1
varrho_bad = varrho*(zh==0)+(zh==1);
% If zh=0 (bad health), there are medical expenses
med_j = (zh==0)*medspend_j*e;
% Out-of-pocket medical expenses
op_j = (1-sub)*med_j;

% calculate gross labor earnings
wage = w*kappa_j*zn*theta_i*varrho_bad*n;

b_tr   = max(0,c_floor + op_j - (1+r)*a - wage - pen_j);
income = wage + pen_j + b_tr - op_j;

% Can deduct medical expenses above a fraction of income:
medical_deduc  = max(0,op_j-deduc*(wage+r*a));

% Taxable income, income taxes and social security taxes
taxable_income = wage+r*a-medical_deduc;
taxes          = fun_tax(taxable_income,tau0,tau1);
taxes_ss       = fun_tax_ss(wage,tau_ss,cap_ss);

% Consumption
cons   = (1+r)*a + income - taxes - taxes_ss - aprime; 

end %end function
