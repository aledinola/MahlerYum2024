function income=Mod_income(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,...
    varrho,w,r,c_floor,sub)
% INPUTS
%  n:       Labor supply (full-time equivalents),d1
%  f:       Health effort, d2
%  aprime:  Next-period assets
%  a:       Current assets
%  zh:      Health shock (0=bad, 1=good), semiz 
%  zn:      Labor productivity shock, z
%  e:       Transitory iid medical shock
%  theta_i: Permanent type (education)
% OUTPUT
%  income:  Income defined as labor earnings plus pensions plus transfers 
%           minus out-of-pocket medical expenses (before taxes)

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

end %end function
