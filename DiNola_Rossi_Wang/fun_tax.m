function tax = fun_tax(y,tau0,tau1)
% INPUTS
%  y:    Taxable income
%  tau0: Level parameter in tax function
%  tau1: Progressivity parameter in tax function
% OUTPUT
%  tax: Income Tax liability

y = max(y,0);

tax = y-tau0*y.^(1-tau1);

% Uncomment the following line if you want to rule out negative taxes
%tax = max(tax,0);

end %end function