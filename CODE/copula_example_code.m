% Example code for some copula functions
%
%  Andrew Patton
%
%  27 February 2006

load 'c:\core\teaching\lse_teaching\lecture_notes\ibm_ccola_rets.txt' -ascii;
ibm = ibm_ccola_rets(:,1);
ccola = ibm_ccola_rets(:,2);


% exceedence correlations
inc = 0.025;
qq = (0.1:inc:0.9)';
qq2 = [(0.1:inc:0.5)';(0.5:inc:0.9)'];

temp1 = ang_chen1(ibm,ccola,qq);
figure(1),plot(qq2,temp1,'o-')


% quantile dependence
temp2 = quantiledep(ibm,ccola,qq);
figure(2),plot(qq,temp2,'o-');


% obtaining Unif(0,1) data
% NOTE: for this example I will just use the empirical cdf to transform the
% data, but in practice this is the step where a model for the conditional
% (marginal) densities would be used: conditional means, variances and
% distributions.

u = empiricalcdf(ibm);
v = empiricalcdf(ccola);
T = length(u);
% dropping the first few so we have an even 2500 observations 
u = u(end-2499:end);
v = v(end-2499:end);
T = length(u)

% estimating some copula models

options = optimset('Display','iter','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6);


% 1. Normal Copula
kappa1 = corrcoef12(norminv(u),norminv(v));
LL1 = NormalCopula_CL(kappa1,[u,v]);	   

% 2. Clayton's copula
lower = 0.0001;
theta0 = 1;
[ kappa2 LL2] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);

% 3. Rotated Clayton copula (with tail dep in upper tail instead of lower)
lower = 0.0001;
theta0 = 1;
[ kappa3 LL3] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,1-[u,v]);

% 4. Plackett copula
lower = 0.0001;
theta0 = 1;
[ kappa4 LL4] = fmincon('plackettCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);
% LL5 = -3.2721

% 5. Frank copula
theta0 = 1;
[ kappa5 LL5] = fmincon('frankCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);

% 6. Gumbel copula
lower = 1.1;
theta0 = 2;
[ kappa6 LL6] = fmincon('gumbelCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);

% 7. Rotated Gumbel copula
lower = 1.1;
theta0 = 2;
[ kappa7 LL7] = fmincon('gumbelCL',theta0,[],[],[],[],lower,[],[],options,1-[u,v]);

% 8. Student's t copula
lower = [-0.9 , 2.1 ];
upper = [ 0.9 , 100 ];
theta0 = [kappa1;10];
[ kappa8 LL8] = fmincon('tcopulaCL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);

% 9. Symmetrised Joe-Clayton copula
lower = [0 , 0 ];
upper = [ 1 , 1];
theta0 = [0.25;0.25];
[ kappa9 LL9] = fmincon('sym_jc_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);

LL = [LL1;LL2;LL3;LL4;LL5;LL6;LL7;LL8;LL9];
[(1:length(LL))',LL]
sortrows([(1:length(LL))',LL],2)

% optimal copula (in terms of log-likelihood) is one with lowest likelihood
% (since we minimise the *negative* LL, rather than maximise the positive LL) 
opt_copula = find(LL==min(LL))
% for these assets it is copula 8, the t-copula, 
% followed by copula 9, the SJC copula
% the worst is copula 3, the rotated Clayton

% tail dependence implied by each of these copulas
tauLU = nines(9,2);
tauLU(1,:) = [0,0];                 % Normal copula has zero tail dependence
tauLU(2,:) = [2^(-1/kappa2),0];     % Clayton copula has zero upper tail dependence
tauLU(3,:) = [0,2^(-1/kappa3)];     % Rotated Clayton copula has zero lower tail dependence
tauLU(4,:) = [0,0];                 % Plackett copula has zero tail dependence
tauLU(5,:) = [0,0];                 % Frank copula has zero tail dependence
tauLU(6,:) = [0,2-2^(1/kappa6)];    % Gumbel copula has zero lower tail dependence
tauLU(7,:) = [2-2^(1/kappa7),0];    % Rotated Gumbel copula has zero upper tail dependence
tauLU(8,:) = ones(1,2)*2*tdis_cdf(-sqrt((kappa8(2)+1)*(1-kappa8(1))/(1+kappa8(1))),kappa8(2)+1);  % Student's t copula has symmetric tail dependence
tauLU(9,:) = kappa9([2,1])';               % SJC copula parameters are the tail dependence coefficients, but in reverse order.
tauLU
% the tail dependence values are reasonably similar, when they are allowed
% to be non-zero
sortrows([(1:9)',LL,tauLU],2)
% the 3 best fitting copulas all allow for non-zero lower tail
% dependence. however, Clayton's copula does poorly even with lower tail
% dependence, suggesting that it is just a poor parameterisation for these
% two stocks.


% Now taking a look at a couple of time-varying copulas

% 10. Time-varying normal Copula
lower = -5*ones(3,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
upper = 5*ones(3,1);
theta0 = [log((1+kappa1)/(1-kappa1));0;0];
[ kappa10 LL10] = fmincon('bivnorm_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa1);
[LL10, rho10] = bivnorm_tvp1_CL(kappa10,[u,v],kappa1);
figure(10),plot((1:T)',rho10,(1:T)',kappa1*ones(T,1),'r--'),legend('time-varying','constant'),title('Normal copula');
% looks nice

% 11. Time-varying rotated Gumbel copula

lower = -5*ones(3,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
upper =  5*ones(3,1);
theta0 = [sqrt(kappa7-1);0;0];
[ kappa11 LL11] = fmincon('Gumbel_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,[1-u,1-v],kappa7);
[LL11, rho11] = Gumbel_tvp1_CL(kappa11,[1-u,1-v],kappa7);
figure(11),plot((1:T)',rho11,(1:T)',kappa7*ones(T,1),'r--'),legend('time-varying','constant'),title('Rotated Gumbel copula');
% not so nice: variation in this parameter looks like just noise.

% 12. Time-varying SJC copula
lower = -25*ones(6,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
upper =  25*ones(6,1);
theta0 = [log(kappa9(1)/(1-kappa9(1)));0;0;log(kappa9(2)/(1-kappa9(2)));0;0];
[ kappa12 LL12] = fmincon('sym_jc_tvp_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa9);
[ LL12 tauU12 tauL12] = sym_jc_tvp_CL(kappa12,[u,v],kappa9);
figure(12),subplot(2,1,1),plot((1:T)',tauL12,(1:T)',kappa9(2)*ones(T,1),'r--'),legend('time-varying','constant'),title('SJC copula - lower tail'),axis([0,T,0,0.8]);
subplot(2,1,2),plot((1:T)',tauU12 ,(1:T)',kappa9(1)*ones(T,1),'r--'),legend('time-varying','constant'),title('SJC copula - upper tail'),axis([0,T,0,0.8]);
% movement in upper tail dependence seems very noisy, whereas some of the movement in lower tail dependence appears informative.

LL = [LL1;LL2;LL3;LL4;LL5;LL6;LL7;LL8;LL9;LL10;LL11;LL12];
[(1:length(LL))',LL]
sortrows([(1:length(LL))',LL],2)
% new rankings:
% 1 is time-varying SJC copula
% 2 is time-varying rotated Gumbel
% 3 is constant Student's t
% 4 is time-varying Normal

params = [ones(7,1);2;2;3;3;6];  % number of parameters in each model
AIC = 2*LL + 2/T*params;
BIC = 2*LL + log(T)/T*params;
[(1:length(LL))',LL,AIC,BIC]
sortrows([(1:length(LL))',LL,AIC,BIC],2)
sortrows([(1:length(LL))',LL,AIC,BIC],3)
sortrows([(1:length(LL))',LL,AIC,BIC],4)
% rankings by AIC and BIC are the same as by log-likelihood (T is so large
% that k=1 vs k=6 does not impose a very large penalty)

