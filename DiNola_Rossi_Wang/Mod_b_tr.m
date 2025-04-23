function b_tr=Mod_b_tr(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,medspend_j,...
    varrho,w,r,c_floor,sub)
% PURPOSE
%  Mod_b_tr computes welfare transfers b_tr
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
%  b_tr:    Welfare transfer to make sure consumption>=c_floor

% If zh=0 (bad health), there is an income penalty. Note that varrho<1
varrho_bad = varrho*(zh==0)+(zh==1);
% If zh=0 (bad health), there are medical expenses
med_j = (zh==0)*medspend_j*e;
% Out-of-pocket medical expenses
op_j = (1-sub)*med_j;

% calculate gross labor earnings
wage = w*kappa_j*zn*theta_i*varrho_bad*n;

b_tr = max(c_floor + op_j - (1+r)*a - wage - pen_j, 0);


end %end function "Mod_b_tr"
