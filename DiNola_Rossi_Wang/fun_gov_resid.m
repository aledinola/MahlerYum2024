function resid = fun_gov_resid(tau0,Params,Gr,flag,vfoptions,simoptions,FnsToEval)
% fun_gov_resid computes the difference between the G in the governmnet
% budget constraint and its target value from the benchmark economy

Params.tau0 = tau0;

govnet = fun_govnet(Params,Gr,flag,vfoptions,simoptions,FnsToEval);

resid = govnet-Params.govnet_bench;

disp('*******************************************************************')
disp('Level parameter in income tax function:')
disp(tau0)
disp('Government budget residual:')
disp(resid)
disp('*******************************************************************')

end %end function "fun_gov_resid"