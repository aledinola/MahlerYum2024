function labincome=Mod_labincome(n,f,aprime,a,zh,zn,e,theta_i,kappa_j,pen_j,varrho,w)
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
%  labincome: defined as labor earnings plus pensions (before taxes)

% If zh=0 (bad health), there is an income penalty. Note that varrho<1
varrho_bad = varrho*(zh==0)+(zh==1);
labincome  = w*kappa_j*zn*theta_i*varrho_bad*n + pen_j;

end %end function
