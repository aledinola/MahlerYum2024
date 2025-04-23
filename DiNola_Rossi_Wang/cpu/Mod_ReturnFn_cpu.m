function F=Mod_ReturnFn_cpu(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,...
    varrho,w,r,c_floor,deduc,tau0,tau1,sub,tau_ss,cap_ss,sigma,delta,nu,gamma,iota,psi)
% INPUTS
%  n:       Labor supply (full-time equivalents), d1
%  f:       Health effort, d2
%  aprime:  Next-period assets, a'
%  a:       Current assets, a
%  zh:      Health shock (0=sick, 1=healthy), semiz
%  zn:      Labor productivity shock, z
%  e:       Transitory iid medical shock, e
%  theta_i: Permanent type (education), PT
% OUTPUT
%  F:       Current-period payoff F(a',a,z,e)+beta*E[V(a',z',e')]

% If zh=0 (bad health), there is an income penalty. Note that varrho<1
varrho_bad = varrho*(zh==0)+(zh==1);
% If zh=0 (bad health), there are medical expenses
med_j = (zh==0)*medspend_j*e;
% Out-of-pocket medical expenses
op_j = (1-sub)*med_j;

% calculate gross labor earnings
wage = w*kappa_j*zn*theta_i*varrho_bad*n;

b_tr   = max(c_floor + op_j - (1+r)*a - wage - pen_j, 0);
income = wage + pen_j - op_j + b_tr;

% Can deduct medical expenses above a fraction of income:
medical_deduc  = max(0,op_j-deduc*(wage+r*a));

% Taxable income, income taxes and social security taxes
taxable_income = wage+r*a-medical_deduc;
taxes          = fun_tax(taxable_income,tau0,tau1);
taxes_ss       = fun_tax_ss(wage,tau_ss,cap_ss);

% Consumption
cons   = (1+r)*a + income - taxes - taxes_ss - aprime; 

% Preference shifter: If zh=0 (bad health), shifter=1-delta, otherwise
% shifter=1
shifter = 1-delta*(zh==0);

F = -inf(size(cons));
% Note preference shifter
indx = cons>0;
effort_disutil = iota*f^(1+1/psi)/(1+1/psi);
F(indx) = (cons(indx).^(1-sigma))/(1-sigma)-shifter*nu*n^(1+1/gamma)/(1+1/gamma)-effort_disutil; 

end %end function