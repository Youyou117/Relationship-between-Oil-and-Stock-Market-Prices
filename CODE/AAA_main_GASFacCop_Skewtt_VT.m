% This code illustrates the estimation of a skew t-t factor copula with GAS 
% dynamics on a sample of 10 variables using "variance targetting" method. 
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

%%% Nodes and weight for Gauss-Legandre quadrature 
Npt_Quad = 50;               
[Wa, Wb] = GLNodeWt(Npt_Quad);
GLweight = [Wa, Wb];

%%% Factor loadings at t = 1
cr_65 = corr(copula_data(1:65,:),'type','spearman');
cr_full = corr(copula_data,'type','spearman');

lam_ini   = rhobar2betabar(cr_65) ;
lam_bar = rhobar2betabar(cr_full) ;

%%% Bound for parameters [alpha, beta, nuinv_z, nuinv_eps, psi_z] 
options = optimset('Display','iter','MaxIter',1000,'TolFun',1e-6,'TolX',1e-6);
lower     =  [ 0.0001 ; 0.01;      0.008;  0.008; -0.6 ];
upper     =  [ 0.2    ; 0.999999;  0.42;   0.42;   0.6 ]; 

%%% Estimating skew t - t factor copula with GAS recursion
tic

theta_ini =  [  0.05  ; 0.95 ;  0.1;  0.1; 0.1];
[theta_NMLE, fobj, exitflag, output1]= fminsearchbnd('LL_GASFacCop_Skewtt_VT', theta_ini, lower, upper, options, copula_data, GLweight, lam_bar, lam_ini);
[obj, LL, lam, log_lam, s] = LL_GASFacCop_Skewtt_VT(theta_NMLE, copula_data, GLweight, lam_bar, lam_ini);

esti_time = toc

save([output_path,'results_VT_from_GASFactorToolbox2.mat']);

