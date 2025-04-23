function taxable_income=Mod_taxable_income(n,f,aprime,a,zh,zn,e,theta_i,...
    kappa_j,pen_j,medspend_j,varrho,w,r,c_floor,deduc,sub)
% INPUTS
%  n:       Labor supply (full-time equivalents),d1
%  f:       Health effort, d2
%  aprime:  Next-period assets, a'
%  a:       Current assets
%  zh:      Health shock (0=sick, 1=healthy)
%  zn:      Labor productivity shock
%  e:       Transitory iid medical shock
%  theta_i: Permanent type (education)
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
%income = wage + pen_j - op_j + b_tr;

% Can deduct medical expenses above a fraction of income:
medical_deduc  = max(0,op_j-deduc*(wage+r*a));

% Taxable income, income taxes and social security taxes
taxable_income = wage+r*a-medical_deduc;

end %end function "Mod_taxable_income"
