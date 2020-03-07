% This code illustrates the estimation of a skew t-t factor copula with GAS 
% dynamics on a sample of 10 variables. 
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 April 2015
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence 
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University. 


input_path = 'D:\Core\Matlab_programs\Factor_copulas\GAS_factor_copula_toolbox\';  % where the copula data was saved
output_path = 'D:\Core\Matlab_programs\Factor_copulas\GAS_factor_copula_toolbox\'; % where the results MAT file should be saved

%%% Load Data;  (T by N) matrix. This data should be copula data, i.e. all elements are in [0,1]
load([input_path,'sample_data.mat'],'copula_data');
[TT, NN] = size(copula_data);

%%% Group code; each firm has one group code
group_code = [1 1 2 2 3 3 4 4 5 5] ;  
Ngroup = max(group_code);

%%% Nodes and weight for Gauss-Legandre quadrature 
Npt_Quad = 50;               
[Wa, Wb] = GLNodeWt(Npt_Quad);
GLweight = [Wa, Wb];

%%% Factor loadings at t = 1
corr_ini = corr(copula_data(1:65,:));
rho_mean = nan(Ngroup,1);

for i = 1:Ngroup
    inx = (group_code == i );
    rho_mean(i,1) = mean(nonzeros(triu(corr_ini(inx,inx),1)));
end

loading_ini = sqrt(abs(rho_mean'./(1-rho_mean')));

%%% Bound for parameters [omega's, alpha, beta, nuinv_z, nuinv_eps, psi_z] 
options = optimset('Display','iter','MaxIter',1000,'TolFun',1e-6,'TolX',1e-6);
lower     =  [ -0.3*ones(Ngroup,1); 0.0001 ; 0.01;      0.008;  0.008; -0.6 ];
upper     =  [  0.3*ones(Ngroup,1); 0.2    ; 0.999999;  0.42;   0.42;   0.6 ]; 


%%% Estimating skew t - t factor copula with GAS recursion
tic

theta_ini =  [  zeros(Ngroup,1) ;0.05  ; 0.95 ;  0.1;  0.1; 0.1];
[theta_NMLE, fobj, exitflag, output1]= fminsearchbnd('LL_GASFacCop_Skewtt_NGroup', theta_ini, lower, upper, options, copula_data, GLweight, group_code, loading_ini);
[obj, LL, lam, log_lam, s] = LL_GASFacCop_Skewtt_NGroup(theta_NMLE, copula_data, GLweight, group_code, loading_ini);

esti_time = toc  % takes around 3.95 hours on a 3.4GHz machine with 16GB of memory

save([output_path,'results_from_GASFactorToolbox2.mat']);

