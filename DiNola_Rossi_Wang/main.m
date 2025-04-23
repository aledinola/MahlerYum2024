%% Life-Cycle Model with health shocks
% Files to write model equations: they start with suffix Mod_
% - Mod_ReturnFn
% - Mod_labincome
% - Mod_income
% - Mod_cons
% Files to run the program on the cpu:
% - VFI_cpu, VFI_onestep
% - Mod_ReturnFn_cpu
% - f_StatDist_cpu
% - mu_onestep
% - reshape_VandPolicy
% - fun_params_cpu
clear,clc,close all,format long g
% Laptop
myf = 'C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab';
% Desktop
%myf = 'C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab';
addpath(genpath(myf))
% Add cpu-specific functions
addpath('cpu')
addpath(fullfile('..','tools'))

%% How does VFI Toolkit think about this?
% Ordering of variables: (aprime,a,z(1),z(2),theta,j)
% 1 Endogenous state variable:           a, assets 
% 2 Stochastic exogenous state variable: z, income shock and health shock
% 1 Permanent type (low vs high):        theta
%   Age:                                 j

%% Set parameters, grids and options
[Params,Gr,flag,vfoptions,simoptions,FnsToEval] = set_params_grids();

if flag.do_bench==1
    % Solve model and obtain government budget residual
    %govnet = fun_govnet(Params,Gr,flag,vfoptions,simoptions,FnsToEval);

    [sol,sol_cpu,mom] = solve_model(Params,Gr,flag,vfoptions,simoptions,FnsToEval);
    govnet = mom.tax.govnet;

    fid=fopen(flag.saveG_fname,'wt');  % overwrite
    fprintf(fid,'%.12f \n',govnet);
    fclose(fid);
    
    % Write results to txt
    txt_export(Params,Gr,mom,flag);
    tic
    txt_export_all(Params,Gr,mom,flag);
    toc
else
    % Read from flag.saveG_fname and call it govnet
    Params.govnet_bench = importdata(flag.saveG_fname);

    myfun = @(tau0) fun_gov_resid(tau0,Params,Gr,flag,vfoptions,simoptions,FnsToEval);
    options = optimset('TolX',1e-6);
    [tau0_root,resid] = fzero(myfun,Params.tau0,options);
    Params.tau0 = tau0_root;
    [sol,sol_cpu,mom] = solve_model(Params,Gr,flag,vfoptions,simoptions,FnsToEval);
    % Write results to txt
    txt_export(Params,Gr,mom,flag);
end

%% Plot model inputs and results
if flag.do_plots==1
    plot_model_inputs
    plot_model_results(Params,Gr,mom,flag);
end

%% Checks

%  Unpack from mom
% ave            = mom.ave;
% ave_age        = mom.ave_age;
% ave_health     = mom.ave_health;
% ave_age_health = mom.ave_age_health;

% a1 = ave.assets
% a2 = sum(ave_age.assets.*Params.mewj)
% 
% a1 = ave.assets
% a2 = sum(ave_health.assets.*[ave.frac_badhealth,1-ave.frac_badhealth]')
% 
% a1 = ave_health.assets(1)
% a2 = sum(ave_age_health.assets(1,:).*ave_age.frac_badhealth.*Params.mewj)/sum(ave_age.frac_badhealth.*Params.mewj)
