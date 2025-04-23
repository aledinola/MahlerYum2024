function [P] = fun_params_cpu(Params,theta_c,j_c)
% Given age j_c and permanent type theta_c, unpack relevant fields from Params

if ~isstruct(Params)
    error('Input variable Params must be a struct')
end

P = struct();

%% Params that depend on PT and age
if theta_c==1
    % low type
    P.pen_j = Params.pen_j.low(j_c);
elseif theta_c==2
    % high type
    P.pen_j = Params.pen_j.high(j_c);
else
    error('theta_c out of bounds')
end

%% Params that depend on PT only
if theta_c==1
    P.varrho = Params.varrho.low;
elseif theta_c==2
    P.varrho = Params.varrho.high;
else
    error('theta_c out of bounds')
end

%% Params that depend on age only
P.kappa_j    = Params.kappa_j(j_c);
P.medspend_j = Params.medspend_j(j_c);
P.sj         = Params.sj(j_c);
P.iota_j     = Params.iota_j(j_c);

%% Others
P.beta    = Params.beta;
P.w       = Params.w;
P.r       = Params.r;
P.c_floor = Params.c_floor;
P.deduc   = Params.deduc;
P.tau0    = Params.tau0;
P.tau1    = Params.tau1;
P.sub     = Params.sub;
P.tau_ss  = Params.tau_ss;
P.cap_ss  = Params.cap_ss;
P.sigma   = Params.sigma;
P.delta   = Params.delta;
P.nu      = Params.nu;
P.gamma   = Params.gamma;
P.psi     = Params.psi;

end %end function fun_params_cpu