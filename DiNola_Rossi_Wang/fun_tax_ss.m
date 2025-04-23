function tax = fun_tax_ss(y_ss,tau_ss,cap_ss)
%  Social security tax function
% INPUTS
%  y_ss:   Taxable income
%  tau_ss: Social security tax rate
%  cap_ss: Social security cap
% OUTPUT
%  tax: Income Tax liability

y = max(y_ss,0);

tax = tau_ss*min(y,cap_ss);

end %end function fun_tax_ss