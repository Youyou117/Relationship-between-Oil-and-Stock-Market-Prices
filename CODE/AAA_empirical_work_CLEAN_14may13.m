% Code to replicate some of the empirical results presented in:
%
% Patton, A.J., 2011, Copula Methods for Forecasting Multivariate Time Series, 
%       in Handbook of Economic Forecasting, Volume 2, 
%       G. Elliott and A. Timmermann eds., Elsevier, Oxford. Forthcoming.
%
%  Andrew Patton
%
%  10 May 2011

% NOTES:
%
% 1: This code was written to be run in blocks, not all at once. Each block generates a figure (perhaps several) as well as some tables that
% correspond to those in the chapter. Some additional information may also be generated. I have broken the code into "chapters" (highlighted yellow in
% Matlab's editor), but I suggest running even sub-chapters separately to see what is going on.
%
% 2: Since many p-values and standard errors in this analysis are based on simulation methods, they will vary slightly from the paper depending
% on the random number generator seed. These variations should be slight. 
%
% 3: I have included "warnings" for parts of the code that are very slow. I use tic/toc to keep track of procedures that might be slow, and if
% something takes more than about 30 seconds I flag it with a warning. (Some things are REALLY slow, so you may want to set it running and leave
% overnight.)
%
% 4: I thank Dong Hwan Oh for his help getting this code ready for distribution. (But I take all the blame for remaining bugs.)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING THE PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Matlab_code\ucsd_garch'));                 % adding the "UCSD GARCH" toolbox to the path directory. (This is needed for some parts of the code below.) You can skip this line if you already have this toolbox on your saved directory path
addpath(genpath('C:\Matlab_code\toolbox\Spatial'));            % adding the "Spatial Econometrics" toolbox to the path directory. (This is needed for some parts of the code below.) You can skip this line if you already have this toolbox on your saved directory path
addpath(genpath('C:\Matlab_code\Copula_Handbook_toolbox'));    % adding the path for *this* toolbox. Change this location to wherever you unzipped the toolbox.


data_path = 'C:\Papers\Handbook\Copulas_Forecasting\';         % where the data spreadsheet is saved
save_path = 'C:\Papers\Handbook\Copulas_Forecasting\';         % where you would like any output MAT files saved
save_name = 'output_14may13';                                  % files will be saved as 'save_name_stage_x.mat'  for x=1,2,...,14


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING IN THE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return series are the SP100 and SP600, starting August 17 1995 (SP600 index start date) ending May 20 2011 (Friday before Handbook conference..)
data = xlsread([data_path,'data_sp100_sp600_19950817_20110520.xlsx'],'Sheet1','a2:c3970');
dates1 = data(:,1);  % 19900103 to 20110520
min(dates1)  % 19950817
max(dates1)  % 20110520
rets1 = data(:,2:end);
ret_str = {'SP100','SP600'};
T = length(rets1)  % 3969

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC SUMMARY STATISTICS AND FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table1a = nan(5,2);
table1a(1,:) = mean(rets1);
table1a(2,:) = std(rets1);
table1a(3,:) = skewness(rets1);
table1a(4,:) = kurtosis(rets1);
table1a(5,:) = [corrcoef12(rets1),corrcoef12(empiricalCDF(rets1))];

info.fmt = '%10.3f';
info.cnames = strvcat('SP100','SP600');
info.rnames = strvcat('.','Mean','Std dev','Skewness','Kurtosis','Correl (lin/rnk)');
sprintf(['Table 1a: Summary statistics'])
mprint(table1a,info)


% get dates of first trade day each year to use in plots
jandates = 1
datesYMD = YYYYMMDD2cols3(dates1);
years = unique(datesYMD(:,1));
for yy=2:length(years);
    jandates = [jandates;min(find(datesYMD(:,1)==years(yy)))];
end
jandates = [jandates;length(datesYMD)];
datesYMD(jandates)

prices = 100*exp(cumsum(rets1(:,1)/100));
prices(:,2) = 100*exp(cumsum(rets1(:,2)/100));
figure(191),subplot(2,1,1),plot((1:T)',prices(:,1),'b-','LineWidth',2);hold on;
plot((1:T)',prices(:,2),'r-'),legend('S&P 100','S&P 600',2);grid on;
title('Prices of S&P 100 and S&P 600 indices'),hold off;
set(gca,'XTick',jandates(1:2:end));
set(gca,'XTickLabel',datestr(datenum(datesYMD(jandates(1:2:end),1),datesYMD(jandates(1:2:end),2),datesYMD(jandates(1:2:end),3)),12));
subplot(2,1,2),plot([-15,15],[0,0],'k--',[0,0],[-15,15],'k--','LineWidth',2);hold on;
plot(rets1(:,1),rets1(:,2),'bo'),grid on;
xlabel('S&P 100 return'),ylabel('S&P 600 return'),axis([-13.5,13.5,-13.5,13.5]),title('Daily returns on S&P 100 and S&P 600');

% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_1.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_1.mat'';'];
evalin('base',temp_str);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBTAINING MODELS FOR THE CONDITIONAL MEAN AND VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since this is a chapter on *forecasting*, we will not focus on estimating
% the unconditional copula: need to look at the conditional copula. For
% that, we need to take care of the conditional mean and variance.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First: conditional mean
resids = nan(T,2);
mean_order = nan(2,3);
tic;
for ii=1:2
    [theta1,sig21,vcv1,order1,resids1] = ARMAX_opt(rets1(:,ii),5,5,'BIC');  % takes about 6 seconds per variable
    mean_order(ii,1:2) = order1';
    mean_order(ii,3) = 1-sig21/cov(rets1(:,ii));
    resids(:,ii) = [zeros(max(order1),1);resids1];
    [ii,toc]
end
toc  % takes 10 seconds for each series
mean_order 
% so optimal models are AR(2) and AR(0)

% saving the parameters of the optimal model for use in the KS and CvM tests below
temp1 = ols(rets1(3:end,1),[ones(T-2,1),rets1(2:end-1,1),rets1(1:end-2,1)]);
ARparams1 = temp1.beta; 
ARparams2 = mean(rets1(:,2));

LL = [5;10];
LBpvals = nan(2,length(LL));
for ii=1:2;
    for ll=1:length(LL);
        temp1 = mlag(resids(:,ii),LL(ll));
        temp2 = nwest(resids(LL(ll)+1:end,ii),[ones(T-LL(ll),1),temp1(LL(ll)+1:end,:)]);
        LBpvals(ii,ll) = 1-chi2cdf( temp2.beta(2:end)'*inv(temp2.vcv(2:end,2:end))*temp2.beta(2:end) , LL(ll));
    end
    figure(200+ii),sacf(resids(:,ii),20);title([ret_str{ii},': LB test p-values for L=5,10: ',num2str(LBpvals(ii,1),'%6.2f'),', ',num2str(LBpvals(ii,2),'%6.2f')]);
end
% so no evidence of remaining autocorrelation from optimal models for the mean


% testing for significance of cross-variable lags
% sp100 optimal model is AR(2):

temp = nwest(rets1(3:end,1),[ones(T-2,1),rets1(2:end-1,1),rets1(1:end-2,1),mlag(rets1(3:end,2),5,0)]);
chi2stat = temp.beta(3:end)'*inv(temp.vcv(3:end,3:end))*temp.beta(3:end)
chi2pval = 1-chi2cdf(chi2stat,5)  % 0.13

temp = nwest(rets1(3:end,2),[ones(T-2,1),mlag(rets1(3:end,1),5,0)]);
chi2stat = temp.beta(3:end)'*inv(temp.vcv(3:end,3:end))*temp.beta(3:end)
chi2pval = 1-chi2cdf(chi2stat,5)  % 0.34  
% good: no cross-equation lags are signif, at least for order 5. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second: the conditional variance
% I don't have code that is as clean as for ARMA model. I will considerthese models:
% Const_vol, ARCH(1), GARCH(1,1), GJR-GARCH(1,1,1), GARCH(2,2), GJR-GARCH(2,2,2)

vol_LL = nan(2,6,2); % will use Normal likelihood to compare these vol models
hhat_ALL = nan(T,2,6);
tic;
for ii=1:2
    % model 1: constant volatility
    hhat_ALL(:,ii,1) = cov(resids(:,ii))*ones(T,1);
    
    % model 2: ARCH(1)
    [parameters2, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),1,0,0,'GJRGARCH','NORMAL',[],[cov(resids(:,ii))*(1-0.3);0.3]);parameters2
    hhat_ALL(:,ii,2) = hhat1;
    
    % model 3: GARCH(1,1)
    [parameters3, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),1,0,1,'GJRGARCH','NORMAL',[],[cov(resids(:,ii))*0.05;0.05;0.9]);parameters3
    hhat_ALL(:,ii,3) = hhat1;
    
    % model 4: GJR-GARCH(1,1,1)
    [parameters4, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),1,1,1,'GJRGARCH','NORMAL',[],[parameters3(1:2);0;parameters3(3)]);parameters4
    hhat_ALL(:,ii,4) = hhat1;
    
    % saving these parametesr for the KS and CvM tests below, as this model turns out to win according to the BIC
    if ii==1;
        GARCHparams1 = parameters4;
    elseif ii==2;
        GARCHparams2 = parameters4;
    end
    
    % model 5: ARCH(2)
    [parameters5, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),2,0,0,'GJRGARCH','NORMAL',[],[cov(resids(:,ii))*(1-0.7);0.5;0.2]);parameters5
    hhat_ALL(:,ii,5) = hhat1;
    
    % model 6: GARCH(2,2)
    [parameters6, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),2,0,2,'GJRGARCH','NORMAL',[],[parameters3(1:2);0;parameters3(3);0]);parameters6
    hhat_ALL(:,ii,6) = hhat1;
    
    % model 7: GJR-GARCH(2,2,2)
    [parameters7, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),2,2,2,'GJRGARCH','NORMAL',[],[parameters4(1:2);0;parameters4(3);0;parameters4(4);0]);parameters7
    hhat_ALL(:,ii,7) = hhat1;
    
    % now computing the mean log-like of each of these models (doing it here, rather than using output from multigarch so that I am sure that the constant vol model is being compared in the right way)
    vol_params = [1,2,3,4,3,5,7];  % number of params in the above models
    for mm=1:7;
        vol_LL(ii,mm,1) = -1/2*log(2*pi) -1/2*mean(log(hhat_ALL(:,ii,mm))) - 1/2*mean( (resids(:,ii).^2)./hhat_ALL(:,ii,mm) );
        vol_LL(ii,mm,2) = -2*vol_LL(ii,mm,1) + log(T)/T*vol_params(mm);  % BIC
    end
    
    [ii,toc]
end
toc  % 10 seconds for both assets. using sensible starting values makes this run smoothly.
vol_LL
GARCHparams12 = [GARCHparams1,GARCHparams2];

vol_order = nan(2,3);
hhat_opt = nan(T,2);
orders_all = [[0,0,0];[1,0,0];[1,1,0];[1,1,1];[2,0,0];[2,2,0];[2,2,2]];
for ii=1:2
    temp = find(vol_LL(ii,:,2)==min(vol_LL(ii,:,2)));
    vol_order(ii,:) = orders_all(temp,:);
    hhat_opt(:,ii) = hhat_ALL(:,ii,temp);  % conditional variance from the BIC-optimal model
end
vol_order
% both choose (1,1,1). Nice.
stdresids = resids./sqrt(hhat_opt);


% scatter plots of data and std resids
figure(151),plot([-15,15],[0,0],'k--',[0,0],[-15,15],'k--','LineWidth',2);hold on;
plot(rets1(:,1),rets1(:,2),'bo'),grid on;
xlabel('S&P 100 return'),ylabel('S&P 600 return'),axis([-13.5,13.5,-13.5,13.5]),title('Daily returns on S&P 100 and S&P 600');

figure(152),plot([-15,15],[0,0],'k--',[0,0],[-15,15],'k--','LineWidth',2);hold on;
plot(stdresids(:,1),stdresids(:,2),'bo'),grid on;
xlabel('S&P 100 return'),ylabel('S&P 600 return'),axis([-8,8,-8,8]),title('Standardized residuals for S&P 100 and S&P 600');

table1b = nan(3+4,2);
table1b(1:3,1) = ARparams1;
table1b(1,2) = ARparams2;
table1b(4:7,1) = GARCHparams1;
table1b(4:7,2) = GARCHparams2;

info.fmt = '%10.3f';
info.cnames = strvcat('SP100','SP600');
info.rnames = strvcat('.','phi0','phi1','phi2','w','a','d','b')
sprintf(['Table 1b: Conditional mean and variance parameters'])
mprint(table1b,info)


% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_2.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_2.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SKEW T MODELS FOR THE STANDARDIZED RESIDUALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'DiffMaxChange',Inf,'DiffMinChange',0,'Algorithm','active-set');

outSKEWT = nan(2,2);  % params
lower = [2.1, -0.99];
upper = [Inf, 0.99 ];
theta0 = [6;0];
for ii=1:2;
    theta1 = fmincon('skewtdis_LL',theta0,[],[],[],[],lower,upper,[],options,stdresids(:,ii));
    outSKEWT(ii,:) = theta1';
end
outSKEWT

Uedf = empiricalCDF(stdresids);  % prob integral transforms using the empirical cdf
Uskewt = nan(T,2);
for ii=1:2;
    Uskewt(:,ii) = skewtdis_cdf(stdresids(:,ii),outSKEWT(ii,1),outSKEWT(ii,2));
end
Uall = Uedf;
Uall(:,:,2) = Uskewt;    % so Uall contains nonparam U's first, then SKEWt U's
size(Uall)

% Plot of histogram of returns against fitted density, and QQ plot
figure(208);
bins=250;
ee = (-8:0.01:8)';
ii=1;[temp1,temp2] = hist(stdresids(:,ii),bins);
subplot(2,2,1),bar(temp2,temp1/(T*(max(stdresids(:,ii))-min(stdresids(:,ii)))/bins),'c');hold on;plot(ee,skewtdis_pdf(ee,outSKEWT(ii,1),outSKEWT(ii,2)),'r','LineWidth',2),title('S&P 100'),xlabel('x'),ylabel('f(x)'),legend('Data','Fitted skew t',2),axis([min(ee),max(ee),0,0.55]); grid on; hold off;
ii=2;[temp1,temp2] = hist(stdresids(:,ii),bins);
subplot(2,2,2),bar(temp2,temp1/(T*(max(stdresids(:,ii))-min(stdresids(:,ii)))/bins),'c');hold on;plot(ee,skewtdis_pdf(ee,outSKEWT(ii,1),outSKEWT(ii,2)),'r-','LineWidth',2),title('S&P 600'),xlabel('x'),ylabel('f(x)'),axis([min(ee),max(ee),0,0.55]); grid on; hold off;
ii=1;qq = sort(empiricalCDF(stdresids(:,ii)));
subplot(2,2,3),plot(skewtdis_inv(qq,outSKEWT(ii,1),outSKEWT(ii,2)),sort(stdresids(:,ii)),'bo');hold on;plot([-10,10],[-10,10],'r--','LineWidth',2),axis([-6.5,4.5,-6.5,4.5]),xlabel('Model quantile'),ylabel('Empirical quantile');grid on;
ii=2;qq = sort(empiricalCDF(stdresids(:,ii)));
subplot(2,2,4),plot(skewtdis_inv(qq,outSKEWT(ii,1),outSKEWT(ii,2)),sort(stdresids(:,ii)),'bo');hold on;plot([-10,10],[-10,10],'r--','LineWidth',2),axis([-6.5,4.5,-6.5,4.5]),xlabel('Model quantile'),ylabel('Empirical quantile');grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implementing a PEE-adjusted KS and CvM test for these two density models

% test statistics
KSstat1  = max(abs(sort(Uskewt(:,1)) - (1:T)'/T))
CvMstat1 = sum(  (sort(Uskewt(:,1)) - (1:T)'/T).^2 )
KSstat2  = max(abs(sort(Uskewt(:,2)) - (1:T)'/T))
CvMstat2 = sum(  (sort(Uskewt(:,2)) - (1:T)'/T).^2 )
    
% using a simulation method to get critical values that account for parameter estimation error (PEE)
% WARNING: THIS PART IS PRETTY SLOW
reps=100;
tic;out1a = dist_test_PEE(T,reps,ARparams1,GARCHparams1,outSKEWT(1,:)',GARCHparams1);toc  % takes 5.87 seconds for reps=10, 60.1 seconds for reps=100, 10.2 mins for reps=1000
tic;out1b = dist_test_PEE(T,reps,ARparams2,GARCHparams2,outSKEWT(2,:)',GARCHparams2);toc  

reps=1000;
tic;out1a = dist_test_PEE(T,reps,ARparams1,GARCHparams1,outSKEWT(1,:)',GARCHparams1);toc  
tic;out1b = dist_test_PEE(T,reps,ARparams2,GARCHparams2,outSKEWT(2,:)',GARCHparams2);toc  

KSpval1 = mean(out1a(:,1)>=KSstat1)  % 0.17, 0.124, 0.140
CvMpval1 = mean(out1a(:,2)>=CvMstat1)  % 0.13, 0.093, 0.106
KSpval2 = mean(out1b(:,1)>=KSstat2)  % 0.49, 0.479
CvMpval2 = mean(out1b(:,2)>=CvMstat2)  % 0.20, 0.222

table1c = nan(2+2,2);
table1c(1:2,:) = outSKEWT';
table1c(3:4,:) = [[KSpval1;CvMpval1],[KSpval2;CvMpval2]];
info.fmt = '%10.3f';
info.cnames = strvcat('SP100','SP600');
info.rnames = strvcat('.','nu','lambda','KS pval','CvM pval')
sprintf(['Table 1c: Skew t density parameters and tests'])
mprint(table1c,info)



% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_3.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_3.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEPENDENCE SUMMARY STATISTICS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rank correlation
tic;[rhoShat] = rankcorrel(stdresids,0.1,1000);toc  % 0.5 seconds for bootreps=1000
rhoShat(1:3)'  % estimate and 90% confidence interval bounds: 0.781, 0.769, 0.793

% quantile dependence
inc = 0.025;
QQ = (inc:inc:1-inc)';
QQ2 = [(inc:inc:0.5)';(0.5:inc:1-inc)'];
bootreps=1000;
tic;[out1qq,out2qq,out3qq] = quantiledep(stdresids(:,1),stdresids(:,2),QQ,bootreps);toc % takes 2.83 seconds for bootreps=1000
figure(131),plot(QQ2,out1qq(:,1),'bo-',QQ2,out1qq(:,1)+1.645*out1qq(:,2),'r--',QQ2,out1qq(:,1)-1.645*out1qq(:,2),'r--','LineWidth',2),...
    legend('quantile dep','90% CI'),xlabel('Quantile (q)'),ylabel('Quantile dep'),title('Quantile dependence for SP100 and SP600 std resids'),...
    xlabel('quantile (q)');grid on;
temp = eye(length(QQ2)/2);
R = [-temp,temp(:,end:-1:1)];  % gives me upper tail minus lower tail, from smallest quantile to 0.5
temp2 = sqrt(diag(R*out2qq*R'));
figure(141),plot(QQ2(1:end/2),R*out1qq(:,1),'bo-',QQ2(1:end/2),R*out1qq(:,1)+1.645*temp2,'r--',QQ2(1:end/2),R*out1qq(:,1)-1.645*temp2,'r--','LineWidth',2);hold on;
plot([0,0.5],[0,0],'k--'),...
    title('Difference in upper and lower quantile dependence'),ylabel('Upper minus lower'),xlabel('quantile (q)'),legend('Estimate','90% CI',4);grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimating tail dependence

%%%%%%%%%%%%%%% NONPARAMETRIC ESTIMATION OF TAIL DEP, BASED ON FRAHM, JUNKER AND SCHMIDT (2005, IME)

bootreps = 1000;
alpha=0.10;
tic;[tauLUhjs, ~, tauBOOThjs] = nonparam_tail_dep(Uedf,[],bootreps,alpha,1);toc  % takes 32 seconds for bootreps=1000
tauLUhjs(:,1)'  %   0.4112    0.1122    0.6636
tauLUhjs(:,2)'  %   0.2295    0.0221    0.5293

% testing for difference in upper and lower tail dep
[tauLUhjs(1,2)-tauLUhjs(1,1),quantile(tauBOOThjs(:,2)-tauBOOThjs(:,1),[0.05,0.95])]  %      -0.1816   -0.5034    0.2530
(tauLUhjs(1,2)-tauLUhjs(1,1))/std(tauBOOThjs(:,2)-tauBOOThjs(:,1))  % -0.7954
mean((tauLUhjs(1,2)-tauLUhjs(1,1))<(tauBOOThjs(:,2)-tauBOOThjs(:,1)))  % 0.595
% not signif. lower tail has stronger tail dep, but cannot conclude it is signif greater than upper tail

% now creating new figure using this estimate of tail dep 
figure(131011),subplot(2,1,1),plot(QQ2,out1qq(:,1),'bo-',QQ2,out1qq(:,1)+1.645*out1qq(:,2),'r--',...
    [0,0],tauLUhjs(1,1)*ones(1,2),'bp',[0,min(QQ2)],[tauLUhjs(1,1),out1qq(1,1)],'b:',...
    [1,1],tauLUhjs(1,2)*ones(1,2),'bp',[max(QQ2),1],[out1qq(end,1),tauLUhjs(1,2)],'b:',...
    [0,min(QQ2)],[tauLUhjs(3,1),out1qq(1,1)+1.645*out1qq(1,2)],'m:',...
    [0,min(QQ2)],[tauLUhjs(2,1),out1qq(1,1)-1.645*out1qq(1,2)],'m:',...
    [max(QQ2),1],[out1qq(end,1)+1.645*out1qq(end,2),tauLUhjs(3,2)],'m:',...
    [max(QQ2),1],[out1qq(end,1)-1.645*out1qq(end,2),tauLUhjs(2,2)],'m:',...
    QQ2,out1qq(:,1)-1.645*out1qq(:,2),'r--','LineWidth',2),...
    legend('quantile dep','90% CI',3),xlabel('Quantile (q)'),ylabel('Quantile dep'),title('Quantile dependence for SP100 and SP600 std resids'),...
    xlabel('quantile (q)');grid on;hold off;
temp = eye(length(QQ2)/2);
R = [-temp,temp(:,end:-1:1)];  % gives me upper tail minus lower tail, from smallest quantile to 0.5
temp1 = R*out1qq(:,1);
temp2 = sqrt(diag(R*out2qq*R'));
subplot(2,1,2),plot(QQ2(1:end/2),temp1,'bo-',QQ2(1:end/2),temp1+1.645*temp2,'r--',...
    [0,0],(tauLUhjs(1,2)-tauLUhjs(1,1))*ones(2,1),'bp',[0,min(QQ2)],[(tauLUhjs(1,2)-tauLUhjs(1,1)),temp1(1)],'b:',...
    [0,min(QQ2)],[quantile(tauBOOThjs(:,2)-tauBOOThjs(:,1),0.95),temp1(1)+1.645*temp2(1)],'m:',...
    [0,min(QQ2)],[quantile(tauBOOThjs(:,2)-tauBOOThjs(:,1),0.05),temp1(1)-1.645*temp2(1)],'m:',...
    QQ2(1:end/2),temp1-1.645*temp2,'r--','LineWidth',2);grid on;hold on;
plot([0,0.5],[0,0],'k--'),...
    title('Difference in upper and lower quantile dependence'),ylabel('Upper minus lower'),xlabel('quantile (q)');hold off;


qstar = 0.025;
bootreps=100;
tic;[tauLhat,tauUhat,thetahatLhat,thetahatUhat,nobs,outboot] = tail_copula_Gumbel(Uedf,qstar,bootreps);toc  % takes 36 seconds for bootreps=100. takes ?? seconds for bootreps=1000
%%%%%%%%%%%%%%%%%%%%%%
% 14may2013: Note that the estimates of tail dependence from the Gumbel copula from this code are different to those reported in the
% published version of the chapter. I revised the code slightly after the chapter went to press, and was not able to update that table.
% The estimates I report in comment format in the following two lines are the correct "updated" estimates.
%%%%%%%%%%%%%%%%%%%%%%
[tauLhat,quantile(outboot(:,1),[0.05,0.95])]  %      0.3839    0.3045    0.4608
[tauUhat,quantile(outboot(:,2),[0.05,0.95])]  %      0.2098    0.1527    0.2712
mean( abs(tauUhat-tauLhat)>=abs(outboot(:,2)-outboot(:,1)) )  % 0.410 => fail to reject null that tail dep is equal


bootreps=100;
tic;[tauLhat,tauUhat,thetahatLhat,thetahatUhat,nobs,outboot] = tail_copula_t(Uedf,qstar,bootreps);toc  % takes 1.24 hours for bootreps=100. SLOW!
[tauLhat,quantile(outboot(:,1),[0.05,0.95])]  %          0.2657    0.2213    0.3489
[tauUhat,quantile(outboot(:,2),[0.05,0.95])]  %        0.1491    0.0813    0.1698
mean( abs(tauUhat-tauLhat)>=abs(outboot(:,2)-outboot(:,1)) )  % 0.230 => fail to reject null that tail dep is equal


%%%% tail dependence based on parametric copulas
QQ = (0.02:0.02:0.1)';
QQ = (0.005:0.005:0.2)';
tauSKEWT = nan(length(QQ),6);
tic;
for qq=1:length(QQ);
    [tauL,tauU,thetahatL,thetahatU,nobs] = tail_copula_Gumbel(Uskewt,QQ(qq));
    tauSKEWT(qq,:) = [tauL,tauU,thetahatL,thetahatU,nobs];
    [qq,toc]
end
toc  % takes 14 seconds for length(QQ)=40
[QQ,tauSKEWT]
figure(123423),plot(QQ,tauSKEWT(:,1:2),'o-'),legend('lower tail','upper tail',2),xlabel('quantile (q)'),ylabel('tail dependence'),title('Estimated tail dependence using Gumbel tail copula');
% clear drift in estimate as we use more obs from the centre of the dist'n. similar to hill estimator. flat spot not immediately obvious, but appears
% to be around 0.015 to 0.03 for lower tail, and around 0.025 to 0.035 for upper tail. Let's use 0.025 as that is in both ranges

qstar=0.025;  
tauLhat = tauSKEWT( (QQ==qstar), 1)  %  0.38948
tauUhat = tauSKEWT( (QQ==qstar), 2)  %  0.26853
thetaLhat = tauSKEWT( (QQ==qstar), 3) %  1.4545
thetaUhat = tauSKEWT( (QQ==qstar), 4)  %  1.2626

bootreps=100;
tic;[tauLhat,tauUhat,thetahatLhat,thetahatUhat,nobs,outboot] = tail_copula_Gumbel(Uskewt,qstar,bootreps);toc  % takes 36 seconds for bootreps=100. takes ?? seconds for bootreps=1000
[tauLhat,quantile(outboot(:,1),[0.05,0.95])]  %     0.3895    0.3213    0.4568
[tauUhat,quantile(outboot(:,2),[0.05,0.95])]  %     0.2685    0.1853    0.3543


tic;[tauLhat,tauUhat,thetahatLhat,thetahatUhat,nobs,outboot] = tail_copula_Gumbel(Uedf,qstar,bootreps);toc  % takes 36 seconds for bootreps=100. takes ?? seconds for bootreps=1000
[tauLhat,quantile(outboot(:,1),[0.05,0.95])]  %     0.3839    0.3066    0.4617
[tauUhat,quantile(outboot(:,2),[0.05,0.95])]  %     0.2098    0.1414    0.2743


figure(13101),subplot(2,1,1),plot(QQ2,out1qq(:,1),'bo-',QQ2,out1qq(:,1)+1.645*out1qq(:,2),'r--',...
    [0,0],tauLhat*ones(1,2),'bp',[0,min(QQ2)],[tauLhat,out1qq(1,1)],'b:',...
    [1,1],tauUhat*ones(1,2),'bp',[max(QQ2),1],[out1qq(end,1),tauUhat],'b:',...
    [0,min(QQ2)],[quantile(outboot(:,1),0.95),out1qq(1,1)+1.645*out1qq(1,2)],'m:',...
    [0,min(QQ2)],[quantile(outboot(:,1),0.05),out1qq(1,1)-1.645*out1qq(1,2)],'m:',...
    [max(QQ2),1],[out1qq(end,1)+1.645*out1qq(end,2),quantile(outboot(:,2),0.95)],'m:',...
    [max(QQ2),1],[out1qq(end,1)-1.645*out1qq(end,2),quantile(outboot(:,2),0.05)],'m:',...
    QQ2,out1qq(:,1)-1.645*out1qq(:,2),'r--','LineWidth',2),...
    legend('quantile dep','90% CI',3),xlabel('Quantile (q)'),ylabel('Quantile dep'),title('Quantile dependence for SP100 and SP600 std resids'),...
    xlabel('quantile (q)');grid on;
temp = eye(length(QQ2)/2);
R = [-temp,temp(:,end:-1:1)];  % gives me upper tail minus lower tail, from smallest quantile to 0.5
temp1 = R*out1qq(:,1);
temp2 = sqrt(diag(R*out2qq*R'));
subplot(2,1,2),plot(QQ2(1:end/2),temp1,'bo-',QQ2(1:end/2),temp1+1.645*temp2,'r--',...
    [0,0],(tauUhat-tauLhat)*ones(2,1),'bp',[0,min(QQ2)],[(tauUhat-tauLhat),temp1(1)],'b:',...
    [0,min(QQ2)],[quantile(outboot(:,2)-outboot(:,1),0.95),temp1(1)+1.645*temp2(1)],'m:',...
    [0,min(QQ2)],[quantile(outboot(:,2)-outboot(:,1),0.05),temp1(1)-1.645*temp2(1)],'m:',...
    QQ2(1:end/2),temp1-1.645*temp2,'r--','LineWidth',2);grid on;hold on;
plot([0,0.5],[0,0],'k--'),...
    title('Difference in upper and lower quantile dependence'),ylabel('Upper minus lower'),xlabel('quantile (q)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testing for asymmetric dependence 

QQ3 = [0.025;0.05;0.10];
QQ3 = [QQ3;sort(1-QQ3)]
[out1a,out2a] = quantiledep(stdresids(:,1),stdresids(:,2),QQ3,bootreps);
temp = eye(length(QQ3)/2);
R = [-temp,temp(:,end:-1:1)];  % gives me upper tail minus lower tail, from smallest quantile to 0.5
chi2stat = (R*out1a(:,1))'*inv(R*out2a*R')*R*out1a(:,1)  % 2.54
chi2pval = 1-chi2cdf(chi2stat,length(QQ3)/2)  % 0.47
% so not significantly different.

% now testing tail dep equality
[tauLhat,quantile(outboot(:,1),[0.05,0.95])]
[tauUhat,quantile(outboot(:,2),[0.05,0.95])]
tstat = (tauUhat-tauLhat)/std(outboot(:,2)-outboot(:,1))  %  -1.9224
2*normcdf(-abs(tstat))   %  0.0546



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rolling window rank correlations

% WARNING: THIS PART IS PRETTY SLOW
rolling_rank_correl = nan(T,3);
window = 60;
bootreps=100;
tic;
for tt=window:T;
    out1a = rankcorrel([stdresids(tt-window+1:tt,1),stdresids(tt-window+1:tt,2)],0.1,bootreps);
    rolling_rank_correl(tt,:) = out1a(1:3)';
end
toc  % takes 137 seconds for bootreps=100

temp124 = ~isnan(rolling_rank_correl(:,1));
figure(182),plot((1:T),rolling_rank_correl(:,1)','k-');hold on;
jbfill(find(temp124)',rolling_rank_correl(temp124,2)',rolling_rank_correl(temp124,3)','r','r',[],0.2);hold on;
plot((1:T)',zeros(T,1),'k--'),title('60-day rolling rank correlation'),legend('Rank correl','90% CI',4);hold off;
title('60-day rolling rank correlation for SP100 and SP600 std resids'),legend('Rank correl','90% CI',4);grid on;hold off;
set(gca,'XTick',jandates(1:2:end));
set(gca,'XTickLabel',datestr(datenum(datesYMD(jandates(1:2:end),1),datesYMD(jandates(1:2:end),2),datesYMD(jandates(1:2:end),3)),12));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testing for time-varying dependence

bootreps=1000;
tic;pval1 = break_test_tv_copula_2(Uedf,0.15,bootreps);toc
tic;pval2 = break_test_tv_copula_2(Uedf,0.50,bootreps);toc
tic;pval3 = break_test_tv_copula_2(Uedf,0.85,bootreps);toc
% takes 0.27 seconds for bootreps=1000

% WARNING: THIS PART IS PRETTY SLOW
bootreps=100
tic;pval4 = break_test_tv_copula_2(Uedf,-0.15,bootreps);toc
bootreps=1000
tic;pval4  = break_test_tv_copula_2(Uedf,-0.15,bootreps);toc
% takes 30 seconds for bootreps=100, about 300 seconds for bootreps=1000


bootreps=1000;
tic;pval5 = AR_test_rank_correl(Uedf,1,bootreps);toc
tic;pval6 = AR_test_rank_correl(Uedf,5,bootreps);toc
tic;pval7 = AR_test_rank_correl(Uedf,10,bootreps);toc
% takes 3.4 seconds for bootreps=1000

table2 = [pval1,pval2,pval3,pval4,pval5,pval6,pval7];
info.fmt = '%10.3f';
info.cnames = strvcat('0.15','0.5','0.85','Anywhere','AR(1)','AR(5)','AR(10)')
info.rnames = strvcat('.','p-value')
sprintf(['Table 2: Tests for time-varying dependence'])
mprint(table2,info)



% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_4.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_4.mat'';'];
evalin('base',temp_str);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATING SOME CONSTANT COPULAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'DiffMaxChange',Inf,'DiffMinChange',0,'Algorithm','active-set');

thetaALL = nan(9,2,2);   % [ copula model] ; [parameters ] ;[EDF or SKEWT marginal dist]
LLALL = nan(9,2);   % [ copula model] ; [EDF or SKEWT marginal dist]
tauLUall = nan(9,2,2);
tic;

for uu=1:2;  % uu=1 uses EDF for marginals, uu=2 uses SKEWT marginal dist
    u = Uall(:,1,uu);
    v = Uall(:,2,uu);
    
    % 1. Normal Copula
    kappa1 = corrcoef12(norminv(u),norminv(v));
    LL1 = NormalCopula_CL(kappa1,[u,v]);
    thetaALL(1,1,uu) = kappa1;
    
    % 2. Clayton's copula
    lower = 0.0001;
    theta0 = 1;
    [ kappa2 LL2] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);
    thetaALL(2,1,uu) = kappa2;
    
    % 3. Rotated Clayton copula (with tail dep in upper tail instead of lower)
    lower = 0.0001;
    theta0 = 1;
    [ kappa3 LL3] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,1-[u,v]);
    thetaALL(3,1,uu) = kappa3;
    
    % 4. Plackett copula
    lower = 0.0001;
    theta0 = 1;
    [ kappa4 LL4] = fmincon('plackettCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);
    thetaALL(4,1,uu) = kappa4;
    
    % 5. Frank copula
    lower = 0.0001;
    upper = 9;  % need to put some loose upper bound here as LL goes to NaN when theta gets too big
    theta0 = 1;
    [ kappa5 LL5] = fmincon('frankCL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
    thetaALL(5,1,uu) = kappa5;
    
    % 6. Gumbel copula
    lower = 1.1;
    upper = 5;
    theta0 = 2;
    [ kappa6 LL6] = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
    thetaALL(6,1,uu) = kappa6;
    
    % 7. Rotated Gumbel copula
    lower = 1.1;
    upper = 5;
    theta0 = 2;
    [ kappa7 LL7] = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,1-[u,v]);
    thetaALL(7,1,uu) = kappa7;
    
    % 8. Student's t copula  (estimating nu_inv rather than nu)
    lower = [-0.9 , 1/100 ];
    upper = [ 0.9 , 1/2.1 ];
    theta0 = [kappa1;10];
    [ kappa8 LL8] = fmincon('tcopulaCL2',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
    thetaALL(8,1:2,uu) = kappa8;
    
    % 9. Symmetrised Joe-Clayton copula
    lower = [0 , 0 ];
    upper = [ 1 , 1];
    theta0 = [0.4;0.4];
    [ kappa9 LL9] = fmincon('sym_jc_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
    thetaALL(9,1:2,uu) = kappa9([2,1]);  % putting tauL before tauU
    
    LLALL(:,uu) = -[LL1;LL2;LL3;LL4;LL5;LL6;LL7;LL8;LL9];
    opt_copula = find(LLALL(:,uu)==max(LLALL(:,uu)))
    
    % tail dependence implied by each of these copulas
    tauLUall(1,:,uu) = [0,0];                 % Normal copula has zero tail dependence
    tauLUall(2,:,uu) = [2^(-1/kappa2),0];     % Clayton copula has zero upper tail dependence
    tauLUall(3,:,uu) = [0,2^(-1/kappa3)];     % Rotated Clayton copula has zero lower tail dependence
    tauLUall(4,:,uu) = [0,0];                 % Plackett copula has zero tail dependence
    tauLUall(5,:,uu) = [0,0];                 % Frank copula has zero tail dependence
    tauLUall(6,:,uu) = [0,2-2^(1/kappa6)];    % Gumbel copula has zero lower tail dependence
    tauLUall(7,:,uu) = [2-2^(1/kappa7),0];    % Rotated Gumbel copula has zero upper tail dependence
    tauLUall(8,:,uu) = ones(1,2)*2*tdis_cdf(-sqrt((kappa8(2)+1)*(1-kappa8(1))/(1+kappa8(1))),kappa8(2)+1);  % Student's t copula has symmetric tail dependence
    tauLUall(9,:,uu) = kappa9([2,1])';               % SJC copula parameters are the tail dependence coefficients, but in reverse order.
end
toc  % takes 108 seconds for 2 pairs and 2 types of U variables, so around 60 seconds per pair
LLALL
tauLUall
LLALLranks = ranks(LLALL)

table3 = nan(9,3+3);
table3(:,[3,6]) = LLALL([1:7,9,8],[2,1]);
table3(:,1:2) = thetaALL([1:7,9,8],1:2,2);
table3(:,4:5) = thetaALL([1:7,9,8],1:2,1);

info.fmt = strvcat('%10.3f','%10.3f','%10.1f','%10.3f','%10.3f','%10.1f');
info.cnames = strvcat('param1','param2','LL','param1','param2','LL');
info.rnames = strvcat('.','Normal','Clayton','Rot Clayton','Plackett','Frank','Gumbel','Rot Gumbel','SJC','Student''s t');
fprintf(['\n\nTable 3: Constant copula model parameter estimates'])
fprintf(['\n                --------- Parametric -------     ------- Semiparametric ----- \n'])
mprint(table3,info)



% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_5.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_5.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARD ERRORS FOR CONSTANT COPULAS - PARAMETRIC CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rather than consider all 9 copulas above, I will focus on just 4 here: 
% Normal (useful benchmark model)
% Clayton (a pretty bad model for this data)
% Rotated Gumbel (a pretty good model for this data)
% Student's t (the best model of the 9 from above)


% storing the parameter estimates and std errors for Normal, Clayton, Rot Gumbel, and Stud t
table4 = nan(5*2,4+4); 
for uu=1:2  % first nonparam, then param margins
    table4(1,-4*(uu-1)+6)  = thetaALL(1,1,uu);
    table4(3,-4*(uu-1)+6)  = thetaALL(2,1,uu);
    table4(5,-4*(uu-1)+6)  = thetaALL(7,1,uu);
    table4(7,-4*(uu-1)+6)  = thetaALL(8,1,uu);
    table4(9,-4*(uu-1)+6)  = thetaALL(8,2,uu);
end
table4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap std errors for two-stage param estimator

% WARNING: THIS PART IS VERY SLOW
bootreps = 100; % this should be at least 100, and better if it is 1000
block_length = 10;
bootdates = stat_bootstrap_function_21(T,bootreps,block_length);
kappa12boot = nan(bootreps,1+1+1+2);
tic;
for bb=1:bootreps
    rets1boot = rets1(bootdates(:,bb),:);
    [~,~,~,resids_1boot] = ARMAX_2(rets1boot(:,1),2,0);
    resids_1boot = [zeros(2,1);resids_1boot];
    [~,~,~,resids_2boot] = ARMAX_2(rets1boot(:,2),0,0);
    [~, ~, ~, ~, hhat1boot] = multigarch_AJP(resids_1boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    [~, ~, ~, ~, hhat2boot] = multigarch_AJP(resids_2boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    stdresids1 = resids_1boot./sqrt(hhat1boot);
    stdresids2 = resids_2boot./sqrt(hhat2boot);
    
    lower = [2.1, -0.99]; upper = [Inf, 0.99 ];
    warning off;
    skewt1boot = fmincon('skewtdis_LL',outSKEWT(1,:)',[],[],[],[],lower,upper,[],options,stdresids1);
    warning off;
    skewt2boot = fmincon('skewtdis_LL',outSKEWT(2,:)',[],[],[],[],lower,upper,[],options,stdresids2);
    Uskewt1 = skewtdis_cdf(stdresids1,skewt1boot(1),skewt1boot(2));
    Uskewt2 = skewtdis_cdf(stdresids2,skewt2boot(1),skewt2boot(2));
    
    % 1. Normal Copula
    kappa1b = corrcoef12(norminv(Uskewt1),norminv(Uskewt2));
    % 2. Clayton's copula -- getting this just in case I want to present it too
    lower = 0.0001;    warning off;
    kappa2b = fmincon('claytonCL',thetaALL(2,1,2),[],[],[],[],lower,[],[],options,[Uskewt1,Uskewt2]);
    % 7. Rotated Gumbel copula
    lower = 1.1;
    upper = 5;
    kappa7b = fmincon('gumbelCL',thetaALL(7,1,2),[],[],[],[],lower,upper,[],options,1-[Uskewt1,Uskewt2]);
    % 8. Student's t copula
    lower = [-0.9 , 0.01 ];
    upper = [ 0.9 , 0.45 ];
    kappa8b = fmincon('tcopulaCL2',thetaALL(8,:,2)',[],[],[],[],lower,upper,[],options,[Uskewt1,Uskewt2]);
    kappa12boot(bb,:) = [kappa1b,kappa2b,kappa7b,kappa8b'];
    
    [bb,toc]  % takes about 5.6 sec per loop
end
toc  % took 92 mins for bootreps=1000. took 19 seconds for bootreps=10
table4(2:2:10,3) = std(kappa12boot);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation std errors for two-stage param estimator

% WARNING: THIS PART IS VERY SLOW
bootreps = 100; % this should be at least 100, and better if it is 1000
kappa12sim = nan(bootreps,1+1+1+2);
GOFsim1 = nan(bootreps,2,4);
GOFsim10 = nan(bootreps,2,4);
tic;
for bb=1:bootreps
    % First simulate from the copula
    U1 = normcdf(mvnrnd(zeros(1,2),[[1,kappa1];[kappa1,1]],T));
    U2 = clayton_rnd(thetaALL(2,1,2),T);
    U3 = 1-Gumbel_rnd(thetaALL(7,1,2),T);
    U4 = tdis_cdf(mvtrnd([[1,thetaALL(8,1,2)];[thetaALL(8,1,2),1]],1/thetaALL(8,2,2),T),1/thetaALL(8,2,2));
    UUU = U1;
    UUU(:,:,2) = U2;
    UUU(:,:,3) = U3;
    UUU(:,:,4) = U4;
    % then obtain the std resids
    EEE = nan(size(UUU));
    for cc=1:4;
        for mm=1:2;
            EEE(:,mm,cc) = skewtdis_inv(UUU(:,mm,cc),outSKEWT(mm,1),outSKEWT(mm,2));
        end
    end
    % next obtain the marginal dynamics
    MMM = nan(size(UUU));  % cond mean
    HHH = nan(size(UUU));  % vol
    YYY = nan(size(UUU));  % simulated raw data
    for cc=1:4
        for mm=1:2;
            HHH(1:2,mm,cc) = GARCHparams12(1,mm)/(1-GARCHparams12(2,mm)-GARCHparams12(3,mm)/2-GARCHparams12(4,mm));   % starting off at unconditional vol
            MMM(1:2,mm,cc) = mean(rets1(:,mm));  % starting off at unconditional mean
            YYY(1:2,mm,cc) = MMM(1:2,mm,cc) + sqrt(HHH(1:2,mm,cc)).*EEE(1:2,mm,cc);
            for tt=3:T
                HHH(tt,mm,cc) = GARCHparams12(1,mm) + GARCHparams12(2,mm)*HHH(tt-1,mm,cc)*(EEE(tt-1,mm,cc)^2) ...
                    + GARCHparams12(3,mm)*HHH(tt-1,mm,cc)*(EEE(tt-1,mm,cc)^2)*(EEE(tt-1,mm,cc)<0) ...
                    + GARCHparams12(4,mm)*HHH(tt-1,mm,cc);
                if mm==1  % then use AR(2) for mean
                    MMM(tt,mm,cc) = ARparams1(1) + ARparams1(2)*YYY(tt-1,mm,cc) + ARparams1(3)*YYY(tt-2,mm,cc);
                else
                    MMM(tt,mm,cc) = ARparams2(1);
                end
                YYY(tt,mm,cc) = MMM(tt,mm,cc) + sqrt(HHH(tt,mm,cc))*EEE(tt,mm,cc);
            end
        end
    end
    
    for cc=1:4;
        % now estimate the models on the simulated data
        [~,~,~,resids_1boot] = ARMAX_2(YYY(:,1,cc),2,0);
        resids_1boot = [zeros(2,1);resids_1boot];
        [~,~,~,resids_2boot] = ARMAX_2(YYY(:,2,cc),0,0);
        [~, ~, ~, ~, hhat1boot] = multigarch_AJP(resids_1boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
        [~, ~, ~, ~, hhat2boot] = multigarch_AJP(resids_2boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
        stdresids1 = resids_1boot./sqrt(hhat1boot);
        stdresids2 = resids_2boot./sqrt(hhat2boot);
        
        lower = [2.1, -0.99]; upper = [Inf, 0.99 ];
        warning off;
        skewt1boot = fmincon('skewtdis_LL',outSKEWT(1,:)',[],[],[],[],lower,upper,[],options,stdresids1);
        warning off;
        skewt2boot = fmincon('skewtdis_LL',outSKEWT(2,:)',[],[],[],[],lower,upper,[],options,stdresids2);
        Uskewt1 = skewtdis_cdf(stdresids1,skewt1boot(1),skewt1boot(2));
        Uskewt2 = skewtdis_cdf(stdresids2,skewt2boot(1),skewt2boot(2));
        
        % later I will want to use a simualtion to obtain critical values for KS and CvM tests of these
        % copula models, so I will combine that step with this one (so that we only run one simulation)
        if cc==1;
            % 1. Normal Copula
            kappa1b = corrcoef12(norminv(Uskewt1),norminv(Uskewt2));
            [KSstat1 CVMstat1] = copula_GOF_stats([Uskewt1,Uskewt2],'NormalCopula_cdf',kappa1b);
            Vhat1 = Uskewt1;  % used in GOF test based on Rosenblatt transforms
            Vhat2 = normcdf(norminv(Uskewt2),kappa1b*norminv(Uskewt1),sqrt(1-kappa1b^2));  % conditional copula of V2 | V1
        elseif cc==2
            % 2. Clayton's copula -- getting this just in case I want to present it too
            lower = 0.0001;    warning off;
            kappa2b = fmincon('claytonCL',thetaALL(2,1,2),[],[],[],[],lower,[],[],options,[Uskewt1,Uskewt2]);
            [KSstat1 CVMstat1] = copula_GOF_stats([Uskewt1,Uskewt2],'clayton_cdf',kappa2b);
            Vhat1 = Uskewt1;
            Vhat2 = ClaytonUgivenV_t(Uskewt2,Uskewt1,0,kappa7b);  % conditional distribution of U2 | U1
        elseif cc==3
            % 7. Rotated Gumbel copula
            lower = 1.1;
            upper = 5;
            kappa7b = fmincon('gumbelCL',thetaALL(7,1,2),[],[],[],[],lower,upper,[],options,1-[Uskewt1,Uskewt2]);
            [KSstat1 CVMstat1] = copula_GOF_stats(1-[Uskewt1,Uskewt2],'gumbel_cdf',kappa7b);
            Vhat1 = 1-Uskewt1;
            Vhat2 = GumbelUgivenV_t(1-Uskewt2,1-Uskewt1,0,kappa7b);  % conditional distribution of U2 | U1
        elseif cc==4
            % 8. Student's t copula
            lower = [-0.9 , 0.01 ];
            upper = [ 0.9 , 0.45 ];
            kappa8b = fmincon('tcopulaCL2',thetaALL(8,:,2)',[],[],[],[],lower,upper,[],options,[Uskewt1,Uskewt2]);
            [KSstat1 CVMstat1] = copula_GOF_stats([Uskewt1,Uskewt2],'tCopula_cdf_new',kappa8b);
            Vhat1 = Uskewt1;
            Vhat2 = mvt_cond_cdf( tinv(Uskewt2,1/kappa8b(2)) , tinv(Uskewt1,1/kappa8b(2)) , zeros(1,2) , [[1,kappa8b(1)];[kappa8b(1),1]] , 1/kappa8b(2));
        end
        GOFsim1(bb,:,cc) = [KSstat1 CVMstat1];
        [KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2],'IndepCop_cdf');  % KS and CVM test on Rosenblatt transforms
        GOFsim10(bb,:,cc) = [KSstat1 CVMstat1];
    end
    kappa12sim(bb,:) = [kappa1b,kappa2b,kappa7b,kappa8b'];
    
    [bb,toc]  % takes about 34 sec per loop
end
toc  % took ?? mins for bootreps=1000. Took ?? seconds for botreps=10

% checking that simulated parameters are roughly around the estimated parameters
[thetaALL(1,1,2),mean(kappa12sim(:,1)),median(kappa12sim(:,1)),min(kappa12sim(:,1)),max(kappa12sim(:,1)),std(kappa12sim(:,1))]
[thetaALL(2,1,2),mean(kappa12sim(:,2)),median(kappa12sim(:,2)),min(kappa12sim(:,2)),max(kappa12sim(:,2)),std(kappa12sim(:,2))]
[thetaALL(7,1,2),mean(kappa12sim(:,3)),median(kappa12sim(:,3)),min(kappa12sim(:,3)),max(kappa12sim(:,3)),std(kappa12sim(:,3))]
[thetaALL(8,:,2)',mean(kappa12sim(:,4:5))',median(kappa12sim(:,4:5))',min(kappa12sim(:,4:5))',max(kappa12sim(:,4:5))',std(kappa12sim(:,4:5))']

table4(2:2:10,4) = std(kappa12sim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive std errors for two-stage param estimator

% Normal copula
uu=2;  % parametric margins 
scoresc = LLgrad_1('NormalCopula_CLa',thetaALL(1,1,uu),Uall(:,1:2,uu));
Hc1 = hessian('NormalCopula_CL',thetaALL(1,1,uu),Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(1,1,uu), sqrt(Vc1naive2)]
table4(2,1)  = sqrt(Vc1naive2);

% Clayton copula
scoresc = LLgrad_1('claytonCLa',thetaALL(2,1,uu),Uall(:,1:2,uu));
Hc1 = hessian('claytonCL',thetaALL(2,1,uu),Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(2,1,uu), sqrt(Vc1naive2) ]
table4(4,1)  = sqrt(Vc1naive2);

% Rot Gumbel copula
scoresc = LLgrad_1('gumbelCLa',thetaALL(7,1,uu),1-Uall(:,1:2,uu));
Hc1 = hessian('gumbelCL',thetaALL(7,1,uu),1-Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(7,1,uu), sqrt(Vc1naive2) ]
table4(6,1)  = sqrt(Vc1naive2);

% Student's t copula
scoresc = LLgrad_1('tcopulaCL2a',thetaALL(8,1:2,uu)',Uall(:,1:2,uu));
Hc1 = hessian('tcopulaCL2',thetaALL(8,1:2,uu)',Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(8,1:2,uu)', sqrt(diag(Vc1naive2))]
table4([8,10],1)  = sqrt(diag(Vc1naive2));
table4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MSMLE std errors for two-stage param estimator

i=1;  j=2;

% first margin stuff
[theta1m,~,~,resids1] = ARMAX_2(rets1(:,i),mean_order(i,1),mean_order(i,2));
resids1 = [zeros(max(mean_order(i,1),mean_order(i,2)),1);resids1];
scores1m = LLgrad_1('ARMA_LLa',theta1m,mean_order(i,1),mean_order(i,2),rets1(:,i));
[theta1v, ~, hessinv1, ~, hhat1, scores1v] = multigarch_AJP(resids1,1,1,1,'GJRGARCH','NORMAL');
hess1v = inv(hessinv1)/T;  % Kevin's code returns the inverse hessian - I convert it back to the hessian here, then divide by T
stdresids1 = resids1./sqrt(hhat1);    
theta1s = outSKEWT(i,:)';
theta1 = [theta1m;theta1v;theta1s];
scores1s = LLgrad_1('skewtdis_LLa',theta1s,stdresids1);
scores1 = [scores1m,scores1v,scores1s];
B1 = newey_west(scores1,floor(T^(1/3)));
H1 = zeros(size(scores1,2),size(scores1,2));
H1m = hessian('ARMA_LL',theta1m,mean_order(i,1),mean_order(i,2),rets1(:,i))/T;
H1v = hess1v;
H1s = hessian_2stage('skewtdis_ARMA_GJRGARCH_LL',theta1,length(theta1s),[],rets1(:,i),mean_order(i,1),mean_order(i,2),1,1,1)/T;
H1(1:length(theta1m),1:length(theta1m)) = H1m;
H1(length(theta1m)+1:length(theta1m)+length(theta1v)  ,length(theta1m)+1:length(theta1m)+length(theta1v)) = H1v;
H1(end-1:end,:) = H1s;
V1 = inv(H1)*B1*(inv(H1)')/T;  % avar[thetahat] = Ainv*B*Ainv, so V[thetahat] ~~ Ainv*B*Ainv/T
[theta1,sqrt(diag(V1)),theta1./sqrt(diag(V1))]

% second margin stuff
[theta2m,~,~,resids2] = ARMAX_2(rets1(:,j),mean_order(j,1),mean_order(j,2));
resids2 = [zeros(max(mean_order(j,1),mean_order(j,2)),1);resids2];
scores2m = LLgrad_1('ARMA_LLa',theta2m,mean_order(j,1),mean_order(j,2),rets1(:,j));
[theta2v, ~, hessinv2, ~, hhat2, scores2v] = multigarch_AJP(resids2,1,1,1,'GJRGARCH','NORMAL');
hess2v = inv(hessinv2)/T;  % Kevin's code returns the inverse hessian - I convert it back to the hessian here, then divide by T
stdresids2 = resids2./sqrt(hhat2);    
theta2s = outSKEWT(j,:)';
theta2 = [theta2m;theta2v;theta2s];
scores2s = LLgrad_1('skewtdis_LLa',theta2s,stdresids2);
scores2 = [scores2m,scores2v,scores2s];
B2 = newey_west(scores2,floor(T^(1/3)));
H2 = zeros(size(scores2,2),size(scores2,2));
H2m = hessian('ARMA_LL',theta2m,mean_order(j,1),mean_order(j,2),rets1(:,j))/T;
H2v = hess2v;
H2s = hessian_2stage('skewtdis_ARMA_GJRGARCH_LL',theta2,length(theta2s),[],rets1(:,j),mean_order(j,1),mean_order(j,2),1,1,1)/T;
H2(1:length(theta2m),1:length(theta2m)) = H2m;
H2(length(theta2m)+1:length(theta2m)+length(theta2v)  ,length(theta2m)+1:length(theta2m)+length(theta2v)) = H2v;
H2(end-1:end,:) = H2s;
V2 = inv(H2)*B2*(inv(H2)')/T;  % avar[thetahat] = Ainv*B*Ainv, so V[thetahat] ~~ Ainv*B*Ainv/T
[theta2,sqrt(diag(V2)),theta2./sqrt(diag(V2))]

% copula stuff:

% Normal copula
uu=2;  % parametric margins 
scoresc = LLgrad_1('NormalCopula_CLa',thetaALL(1,1,uu),Uall(:,1:2,uu));
scoresALL = [scores1,scores2,scoresc];
BB = newey_west(scoresALL,floor(T^(1/3)));
thetac = thetaALL(1,1,uu);
tic;
Hc = hessian_2stage('margin_margin_copula_CL1',...
   [theta1;theta2;thetac],size(thetac,1),[],rets1(:,[i,j]),...
   'skewtdis_ARMA_GJRGARCH_LL','skewtdis_ARMA_GJRGARCH_LL','NormalCopula_CL',...
   size(theta1,1),size(theta2,1),5,5,...
   mean_order(i,1),mean_order(i,2),1,1,1,...
   mean_order(j,1),mean_order(j,2),1,1,1)/T;
toc  % 1.6 seconds
thetaALLmle = [theta1;theta2;thetac];
HH = zeros(length(thetaALLmle),length(thetaALLmle));
HH(1:length(theta1),1:length(theta1)) = H1;
HH(length(theta1) + (1:length(theta2)),length(theta1) + (1:length(theta2))) = H2;
HH(end-length(thetac)+1:end,:) = Hc;
VALL = inv(HH)*BB*(inv(HH)')/T;
[theta1,sqrt(diag(V1)),theta1./sqrt(diag(V1))]
[theta2,sqrt(diag(V2)),theta2./sqrt(diag(V2))]
[thetaALLmle,sqrt(diag(VALL)),thetaALLmle./sqrt(diag(VALL))]
table4(2,2)  = sqrt(VALL(end,end));


% Clayton copula
scoresc = LLgrad_1('claytonCLa',thetaALL(2,1,uu),Uall(:,1:2,uu));
scoresALL = [scores1,scores2,scoresc];
BB = newey_west(scoresALL,floor(T^(1/3)));
thetac = thetaALL(2,1,uu);
tic;
Hc = hessian_2stage('margin_margin_copula_CL1',...
   [theta1;theta2;thetac],size(thetac,1),[],rets1(:,[i,j]),...
   'skewtdis_ARMA_GJRGARCH_LL','skewtdis_ARMA_GJRGARCH_LL','claytonCL',...
   size(theta1,1),size(theta2,1),5,5,...
   mean_order(i,1),mean_order(i,2),1,1,1,...
   mean_order(j,1),mean_order(j,2),1,1,1)/T;
toc  % 1.6 seconds
thetaALLmle = [theta1;theta2;thetac];
HH = zeros(length(thetaALLmle),length(thetaALLmle));
HH(1:length(theta1),1:length(theta1)) = H1;
HH(length(theta1) + (1:length(theta2)),length(theta1) + (1:length(theta2))) = H2;
HH(end-length(thetac)+1:end,:) = Hc;
VALL = inv(HH)*BB*(inv(HH)')/T;
[theta1,sqrt(diag(V1)),theta1./sqrt(diag(V1))]
[theta2,sqrt(diag(V2)),theta2./sqrt(diag(V2))]
[thetaALLmle,sqrt(diag(VALL)),thetaALLmle./sqrt(diag(VALL))]
table4(4,2)  = sqrt(VALL(end,end));

% Rot Gumbel copula
scoresc = LLgrad_1('gumbelCLa',thetaALL(7,1,uu),1-Uall(:,1:2,uu));
scoresALL = [scores1,scores2,scoresc];
BB = newey_west(scoresALL,floor(T^(1/3)));
thetac = thetaALL(7,1,uu);
tic;
Hc = hessian_2stage('margin_margin_copula_CL1',...
   [theta1;theta2;thetac],size(thetac,1),[],rets1(:,[i,j]),...
   'skewtdis_ARMA_GJRGARCH_LL','skewtdis_ARMA_GJRGARCH_LL','rotgumbelCL',...
   size(theta1,1),size(theta2,1),5,5,...
   mean_order(i,1),mean_order(i,2),1,1,1,...
   mean_order(j,1),mean_order(j,2),1,1,1)/T;
toc  % 1.6 seconds
thetaALLmle = [theta1;theta2;thetac];
HH = zeros(length(thetaALLmle),length(thetaALLmle));
HH(1:length(theta1),1:length(theta1)) = H1;
HH(length(theta1) + (1:length(theta2)),length(theta1) + (1:length(theta2))) = H2;
HH(end-length(thetac)+1:end,:) = Hc;
VALL = inv(HH)*BB*(inv(HH)')/T;
[theta1,sqrt(diag(V1)),theta1./sqrt(diag(V1))]
[theta2,sqrt(diag(V2)),theta2./sqrt(diag(V2))]
[thetaALLmle,sqrt(diag(VALL)),thetaALLmle./sqrt(diag(VALL))]
table4(6,2)  = sqrt(VALL(end,end));


% Student's t copula
scoresc = LLgrad_1('tcopulaCL2a',thetaALL(8,1:2,uu)',Uall(:,1:2,uu));
scoresALL = [scores1,scores2,scoresc];
BB = newey_west(scoresALL,floor(T^(1/3)));
thetac = thetaALL(8,1:2,uu)';
tic;
Hc = hessian_2stage('margin_margin_copula_CL1',...
   [theta1;theta2;thetac],size(thetac,1),[],rets1(:,[i,j]),...
   'skewtdis_ARMA_GJRGARCH_LL','skewtdis_ARMA_GJRGARCH_LL','tcopulaCL2',...
   size(theta1,1),size(theta2,1),5,5,...
   mean_order(i,1),mean_order(i,2),1,1,1,...
   mean_order(j,1),mean_order(j,2),1,1,1)/T;
toc  % 1.6 seconds
thetaALLmle = [theta1;theta2;thetac];
HH = zeros(length(thetaALLmle),length(thetaALLmle));
HH(1:length(theta1),1:length(theta1)) = H1;
HH(length(theta1) + (1:length(theta2)),length(theta1) + (1:length(theta2))) = H2;
HH(end-length(thetac)+1:end,:) = Hc;
VALL = inv(HH)*BB*(inv(HH)')/T;
[theta1,sqrt(diag(V1)),theta1./sqrt(diag(V1))]
[theta2,sqrt(diag(V2)),theta2./sqrt(diag(V2))]
[thetaALLmle,sqrt(diag(VALL)),thetaALLmle./sqrt(diag(VALL))]
table4([8,10],2)  = sqrt(diag(VALL(end-1:end,end-1:end)));
table4



% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_6.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_6.mat'';'];
evalin('base',temp_str);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARD ERRORS FOR CONSTANT COPULAS - SEMIPARAMETRIC CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap std errors for two-stage SEMIPARAM estimator

uu=1;  % EDF margins

% WARNING: THIS PART IS PRETTY SLOW
bootreps = 100; % this should be at least 100, and better if it is 1000
bootdates = randint(1,T,T,bootreps);  % iid bootstrap, justified by Chen+Fan (2006) and Remillard (2010)
kappa1boot = nan(bootreps,1+1+1+2);
tic;
for bb=1:bootreps
    Uboot = Uall(bootdates(:,bb),1:2,uu);
    % 1. Normal Copula
    kappa1b = corrcoef12(norminv(Uboot));
    % 2. Clayton's copula -- getting this just in case I want to present it too
    lower = 0.0001;    warning off;
    kappa2b = fmincon('claytonCL',thetaALL(2,1,uu),[],[],[],[],lower,[],[],options,Uboot);
    % 7. Rotated Gumbel copula
    lower = 1.1;
    upper = 5;
    kappa7b = fmincon('gumbelCL',thetaALL(7,1,uu),[],[],[],[],lower,upper,[],options,1-Uboot);
    % 8. Student's t copula
    lower = [-0.9 , 0.01 ];
    upper = [ 0.9 , 0.45 ];
    kappa8b = fmincon('tcopulaCL2',thetaALL(8,:,uu)',[],[],[],[],lower,upper,[],options,Uboot);
    kappa1boot(bb,:) = [kappa1b,kappa2b,kappa7b,kappa8b'];
    [bb,toc]  % takes around 0.7 seconds per replication
end
toc  % takes around 11 mins for bootreps=1000

table4(2,7)  =  nanstd(kappa1boot(:,1));
table4(4,7)  =  nanstd(kappa1boot(:,2));
table4(6,7)  =  nanstd(kappa1boot(:,3));
table4(8,7)  =  nanstd(kappa1boot(:,4));
table4(10,7) =  nanstd(kappa1boot(:,5));

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive std errors for two-stage nonparam estimator

uu=1;  % nonparametric margins 

% Normal copula
scoresc = LLgrad_1('NormalCopula_CLa',thetaALL(1,1,uu),Uall(:,1:2,uu));
Hc1 = hessian('NormalCopula_CL',thetaALL(1,1,uu),Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(1,1,uu), sqrt(Vc1naive2) ]
table4(2,5)  = sqrt(Vc1naive2);

% Clayton copula
scoresc = LLgrad_1('claytonCLa',thetaALL(2,1,uu),Uall(:,1:2,uu));
Hc1 = hessian('claytonCL',thetaALL(2,1,uu),Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(2,1,uu), sqrt(Vc1naive2) ]
table4(4,5)  = sqrt(Vc1naive2);

% Rot Gumbel copula
scoresc = LLgrad_1('gumbelCLa',thetaALL(7,1,uu),1-Uall(:,1:2,uu));
Hc1 = hessian('gumbelCL',thetaALL(7,1,uu),1-Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(7,1,uu), sqrt(Vc1naive2) ]
table4(6,5)  = sqrt(Vc1naive2);

% Student's t copula
scoresc = LLgrad_1('tcopulaCL2a',thetaALL(8,1:2,uu)',Uall(:,1:2,uu));
Hc1 = hessian('tcopulaCL2',thetaALL(8,1:2,uu)',Uall(:,1:2,uu))/T;  % used in naive VCV matrix for copula params
Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
[thetaALL(8,1:2,uu)', sqrt(diag(Vc1naive2))]
table4([8,10],5)  = sqrt(diag(Vc1naive2));
table4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct two-stage std errors for semiparam model 

uu=1;  % nonparam margins
i=1;  j=2;

uS = Uall(:,1,uu);
vS = Uall(:,2,uu);

% Normal copula
thetacS = thetaALL(1,1,uu);
scorescS = LLgrad_1('NormalCopula_CLa',thetacS,[uS,vS]);
Hc1S = hessian('NormalCopula_CL',thetacS,[uS,vS])/T;  % used in naive VCV matrix for copula params. This is Bhat in Chen and Fan (2006)
Vc2S = inv(Hc1S)/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions and from EDF
% now using Chen and Fan (2006) to get correct std errors for this estimate
tic;[sigmaS1,sigmaS2]= copula_deriv2_chen_fan('NormalCopula_CLa',thetacS,[uS,vS]);toc  % takes 2.8 seconds 
VcSa = inv(Hc1S)*sigmaS1*inv(Hc1S)/T;
VcSb = inv(Hc1S)*sigmaS2*inv(Hc1S)/T;
[thetacS,sqrt(diag(Vc2S)),sqrt(diag(VcSa)),sqrt(diag(VcSb))]
table4(2,6) = sqrt(VcSa);

% Clayton copula
thetacS = thetaALL(2,1,uu);
scorescS = LLgrad_1('claytonCLa',thetacS,[uS,vS]);
Hc1S = hessian('claytonCL',thetacS,[uS,vS])/T;  % used in naive VCV matrix for copula params. This is Bhat in Chen and Fan (2006)
Vc2S = inv(Hc1S)/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions and from EDF
% now using Chen and Fan (2006) to get correct std errors for this estimate
tic;[sigmaS1,sigmaS2]= copula_deriv2_chen_fan('claytonCLa',thetacS,[uS,vS]);toc  % takes 2.8 seconds 
VcSa = inv(Hc1S)*sigmaS1*inv(Hc1S)/T;
VcSb = inv(Hc1S)*sigmaS2*inv(Hc1S)/T;
[thetacS,sqrt(diag(Vc2S)),sqrt(diag(VcSa)),sqrt(diag(VcSb))]
table4(4,6) = sqrt(VcSa);

% Rotated Gumbel copula
thetacS = thetaALL(7,1,uu);
scorescS = LLgrad_1('gumbelCLa',thetacS,[uS,vS]);
Hc1S = hessian('gumbelCL',thetacS,[uS,vS])/T;  % used in naive VCV matrix for copula params. This is Bhat in Chen and Fan (2006)
Vc2S = inv(Hc1S)/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions and from EDF
% now using Chen and Fan (2006) to get correct std errors for this estimate
tic;[sigmaS1,sigmaS2]= copula_deriv2_chen_fan('gumbelCLa',thetacS,[uS,vS]);toc  % takes 2.8 seconds 
VcSa = inv(Hc1S)*sigmaS1*inv(Hc1S)/T;
VcSb = inv(Hc1S)*sigmaS2*inv(Hc1S)/T;
[thetacS,sqrt(diag(Vc2S)),sqrt(diag(VcSa)),sqrt(diag(VcSb))]
table4(6,6) = sqrt(VcSa);

% Student's t copula
thetacS = thetaALL(8,1:2,uu)';
scorescS = LLgrad_1('tcopulaCL2a',thetacS,[uS,vS]);
Hc1S = hessian('tcopulaCL2',thetacS,[uS,vS])/T;  % used in naive VCV matrix for copula params. This is Bhat in Chen and Fan (2006)
Vc2S = inv(Hc1S)/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions and from EDF
% now using Chen and Fan (2006) to get correct std errors for this estimate
tic;[sigmaS1,sigmaS2]= copula_deriv2_chen_fan('tcopulaCL2a',thetacS,[uS,vS]);toc  % takes 5.91 seconds 
VcSa = inv(Hc1S)*sigmaS1*inv(Hc1S)/T;
VcSb = inv(Hc1S)*sigmaS2*inv(Hc1S)/T;
[thetacS,sqrt(diag(Vc2S)),sqrt(diag(VcSa)),sqrt(diag(VcSb))]
table4([8,10],6) = sqrt(diag(VcSa));


% simulation-based approach for the constant copulas - semi-parametric
% This is based on Remillard 2010

uu=1;  % nonparametric margins

bootreps = 100;  % should be at least 100
kappaGOFsim3 = nan(bootreps,1+1+1+2);
GOFsim3 = nan(bootreps,2,4);
GOFsim30 = nan(bootreps,2,4);
tic;
for bb=1:bootreps
    % First simulate from the copula
    U1 = normcdf(mvnrnd(zeros(1,2),[[1,thetaALL(1,1,uu)];[thetaALL(1,1,uu),1]],T));
    U2 = clayton_rnd(thetaALL(2,1,uu),T);
    U3 = 1-Gumbel_rnd(thetaALL(7,1,uu),T);
    U4 = tdis_cdf(mvtrnd([[1,thetaALL(8,1,uu)];[thetaALL(8,1,uu),1]],1/thetaALL(8,2,uu),T),1/thetaALL(8,2,uu));
    UUU = U1;
    UUU(:,:,2) = U2;
    UUU(:,:,3) = U3;
    UUU(:,:,4) = U4;
    
    for cc=1:4;
        % then obtain the PITs of the *copula* simulation (seems weird, but this is using the EDF margin estimator)
        Uhat = empiricalCDF( UUU(:,:,cc) );  % using EDF to estimate marginal distributions
        
        if cc==1;
            % 1. Normal Copula
            kappa1b = corrcoef12(norminv(Uhat));
            [KSstat1 CVMstat1] = copula_GOF_stats(Uhat,'NormalCopula_cdf',kappa1b);
            Vhat1 = Uhat(:,1);  % computing the rosenblatt transforms
            Vhat2 = normcdf(norminv(Uhat(:,2)),kappa1b*norminv(Uhat(:,1)),sqrt(1-kappa1b^2));  % conditional copula of V2 | V1
        elseif cc==2
            % 2. Clayton's copula -- getting this just in case I want to present it too
            lower = 0.0001;    warning off;
            theta0 = thetaALL(2,1,uu) + randn/100;
            kappa2b = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,Uhat);
            [KSstat1 CVMstat1] = copula_GOF_stats(Uhat,'clayton_cdf',kappa2b);
            Vhat1 = Uhat(:,1);
            Vhat2 = ClaytonUgivenV_t(Uhat(:,2),Uhat(:,1),0,kappa7b);  % conditional distribution of U2 | U1
        elseif cc==3
            % 7. Rotated Gumbel copula
            lower = 1.1;
            upper = 5;
            theta0 = thetaALL(7,1,uu) + randn/100;
            kappa7b = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,1-Uhat);
            [KSstat1 CVMstat1] = copula_GOF_stats(1-Uhat,'gumbel_cdf',kappa7b);
            Vhat1 = 1-Uhat(:,1);
            Vhat2 = GumbelUgivenV_t(1-Uhat(:,2),1-Uhat(:,1),0,kappa7b);  % conditional distribution of U2 | U1
        elseif cc==4
            % 8. Student's t copula
            lower = [-0.9 , 0.01 ];
            upper = [ 0.9 , 0.45 ];
            theta0 = thetaALL(8,:,uu)' + randn(2,1)/100;
            kappa8b = fmincon('tcopulaCL2',theta0,[],[],[],[],lower,upper,[],options,Uhat);
            [KSstat1 CVMstat1] = copula_GOF_stats(Uhat,'tCopula_cdf_new',kappa8b);
            Vhat1 = Uhat(:,1);
            Vhat2 = mvt_cond_cdf(tinv(Uhat(:,2),1/kappa8b(2)),tinv(Uhat(:,1),1/kappa8b(2)),zeros(1,2),[[1,kappa8b(1)];[kappa8b(1),1]],1/kappa8b(2));
        end
        GOFsim3(bb,:,cc) = [KSstat1 CVMstat1];
        [KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2],'IndepCop_cdf');
        GOFsim30(bb,:,cc) = [KSstat1 CVMstat1];
    end
    kappaGOFsim3(bb,:) = [kappa1b,kappa2b,kappa7b,kappa8b'];
    
    [bb,toc]  % takes about 62 sec per loop, so about 100 mins for 100 replications
end
toc  %
table4([2;4;6;8;10],8) = std(kappaGOFsim3);

table4a = nan(3+3+3+5,4+4); % insert LogL values for each model
table4a(1:2,:) = table4(1:2,:);
table4a(4:5,:) = table4(3:4,:);
table4a(7:8,:) = table4(5:6,:);
table4a(10:13,:) = table4(7:10,:);
table4a([3,6,9,14],[3,6]) = LLALL([1,2,7,8],[2,1]);  % log likelihood values for these models
table4 = table4a;

info.fmt = strvcat('%10.4f');
info.cnames = strvcat('Naive','MSML','Boot','Sim','Naive','MSML','Boot','Sim');
info.rnames = strvcat('.','Normal-rho','Std err','log L','Clayton-kappa','Std err','log L','RotGumbel-kappa','Std err','log L','Stud_t-rho','Std err','Stud_t-nuinv','Std err','log L');
fprintf(['\n\nTable 4: Standard errors on estimated constant copula parameters'])
fprintf(['\n                     -------------- Parametric ------------      ----------- Semiparametric ----------- \n'])
mprint(table4,info)




% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_7.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_7.mat'';'];
evalin('base',temp_str);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATING TIME-VARYING COPULAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I will use the GAS(1,1) specification of Creal et al. (2011) to let the parameter vary through time

table5 = nan(7+9,4+2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-varying rotated Gumbel

% First, it is useful to get the expected Hessian for a range of values of kappa. Only need to do this once, and then can use this in estimation. 
% (See Appendix A of Creal et al. )
% Will do this using simulation

% WARNING: THIS PART IS PRETTY SLOW
reps = 1000;
KK = [(1.1:0.1:2)';(2.25:0.25:4)';5;6];  % values of kappa at which to evaluate the function
KK = [1.5;2;3;4];
KK = [1.01;1.25;(1.5:0.5:5)';6;10];
KK = [1.01;(1.1:0.1:3)';(3.25:0.25:5)';(6:10)'];
HESSgumbel = nan(length(KK),2);
tic;
for ii=1:length(KK);
    Usim = Gumbel_rnd(KK(ii),reps);
    Usim = 1-Usim;  % rotating the data
    scoreSIM = LLgrad_1('gumbelCLa',KK(ii),Usim);
    HESSgumbel(ii,1) = mean(scoreSIM.^2);   
    HESSgumbel(ii,2) = corrcoef12(Usim);  % rank correlation
    [ii,toc]  % takes about 6 seconds per valuke of kappa, for reps=10,000
end
toc  % took 71 seconds for 12 values of kappa and reps=10,000
% took 198 seconds for 34 values of kappa and reps=10,000
figure(3423),plot(KK,HESSgumbel(:,1),'bo-')  % nice and smooth.
figure(34250),plot(KK,HESSgumbel(:,2),'bo-'),xlabel('Gumbel copula parameter'),ylabel('Rank correlation');


% Now estimating the GAS model for the Rotated Gumbel copula
lower = [-3 ; -3 ; 0 ];
upper = [ 3 ; 3 ; 0.999 ];
theta0 = [log(thetaALL(7,1,uu)-1)*(1-0.05-0.93);0.05;0.93];

uu=2;  % skew t margins
tic;kappa8b3 = fmincon('Rotgumbel_GAS_CL',theta0,[],[],[],[],lower,upper,[],options,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);toc  % 69 seconds
kappa8b3'  %      0.0012935     0.040425      0.99609
[CL,kappat,ft] = Rotgumbel_GAS_CL(kappa8b3,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
figure(12311),plot(ft)
figure(12312),plot(kappat)
table5([1;3;5],3) = kappa8b3;
table5(7,3) = -CL;

uu=1;  % EDF margins
tic;kappa7tvs = fmincon('Rotgumbel_GAS_CL',kappa8b3,[],[],[],[],lower,upper,[],options,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);toc  % 50 seconds
kappa7tvs'  %   0.0015    0.0420    0.9955
[CL,kappat,ft] = Rotgumbel_GAS_CL(kappa7tvs,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
table5([1;3;5],5) = kappa7tvs;
table5(7,5) = -CL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-varying Student's t

% First, need to get the expected Hessian for a range of values of [rho,nu]. Only need to do this once, and then can use this in estimation. 
% (See Appendix A of Creal et al. )
% Will do this using simulation
reps = 10000/4;
NN = [2.1;3;4;6;10;20;100];
RR = [-0.99;-0.5;0;0.5;0.99];
RR = [-0.99;-0.75;-0.5;-0.25;-0.1;0;0.1;0.25;0.5;0.75;0.8;0.85;0.9;0.95;0.99];  % adding a few extra nodes for rho around 0.85, which seems to be a common value for this series
HESSstudt = nan(length(RR),length(NN),2,2);
RHOSstudt = nan(length(RR),length(NN),1);
tic;
for rr=1:length(RR);
    for nn=1:length(NN);
        temp = mvtrnd([1,RR(rr);RR(rr),1],NN(nn),reps);
        temp = [temp;-temp];  % imposing symmetry
        temp = [temp;[temp(:,2),temp(:,1)]];  % imposing exchangeability
        U = tdis_cdf(temp(:,1),NN(nn));
        V = tdis_cdf(temp(:,2),NN(nn));
        scoreSIM = LLgrad_1('tcopulaCLa',[RR(rr);NN(nn)],[U,V]);
        HESSstudt(rr,nn,1,1) = mean(scoreSIM(:,1).^2);
        HESSstudt(rr,nn,2,2) = mean(scoreSIM(:,2).^2);
        HESSstudt(rr,nn,1,2) = mean(scoreSIM(:,1).*scoreSIM(:,2));
        HESSstudt(rr,nn,2,1) = HESSstudt(rr,nn,1,2);
        RHOSstudt(rr,nn) = corrcoef12([U,V]);  % rank correlation
        [rr,nn,toc]   % takes about 6 seconds per valuke of kappa, for reps=10,000
    end
end
toc  % took 35 seconds for length(RR)*length(NN)=35 and reps=10,000

% Now estimating the GAS model for the Student's t copula
% WARNING: THIS PART IS VERY SLOW

uu=2;  % skew t margins
lower = [-3 ; -3 ; 0 ; 0.01];
upper = [ 3 ; 3 ; 0.999 ; 0.45];
theta0 = [log( (0.9999+thetaALL(8,1,uu))/(0.9999-thetaALL(8,1,uu)) )*(1-0.05-0.8);0.05;0.8;thetaALL(8,2,uu)];
tic;theta8tv2 = fmincon('tcopula_GAS_CL',theta0,[],[],[],[],lower,upper,[],options,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);toc  % about 30 minuts
theta8tv2'  % 0.019872     0.065285      0.99121     0.088728

tic;[CL,rhot,ft] = tcopula_GAS_CL(theta8tv2,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);toc  % 6.8 secs per function evaluation
table5([8;10;12;14],3) = theta8tv2;
table5(16,3) = -CL;

uu=1;  % EDF margins
tic;theta8tvs = fmincon('tcopula_GAS_CL',theta8tv2,[],[],[],[],lower,upper,[],options,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);toc  % 19 mins
theta8tvs' % 0.0192    0.0603    0.9913    0.0891

tic;[CL,rhot,ft] = tcopula_GAS_CL(theta8tvs,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);toc  % 6.8 secs per functino evaluation
table5([8;10;12;14],5) = theta8tvs;
table5(16,5) = -CL;



% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_8.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_8.mat'';'];
evalin('base',temp_str);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STANDARD ERRORS ON TIME-VARYING COPULAS - PARAMETRIC CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive std errors for param model 

uu=2;  % skew t margins
% Rot Gumbel GAS copula
scoresc = LLgrad_1('Rotgumbel_GAS_CLa',kappa8b3,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
Hc1 =     hessian('Rotgumbel_GAS_CL',  kappa8b3,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)])/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
Bc = newey_west(scoresc,floor(T^(1/3)));
Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
[kappa8b3, sqrt(diag(Vc2)) ]
table5([2;4;6],1) = sqrt(diag(Vc2));

% Student's t GAS copula
% WARNING: THIS PART IS A BIT SLOW
tic;scoresc = LLgrad_1('tcopula_GAS_CLa',theta8tv2,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);toc
tic;Hc1 =      hessian('tcopula_GAS_CL', theta8tv2,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt)/T;toc  % used in naive VCV matrix for copula params
Bc = newey_west(scoresc,floor(T^(1/3)));
Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
[theta8tv2, sqrt(diag(Vc2))]
table5([9;11;13;15],1) = sqrt(diag(Vc2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct MSML std errors for param model 

i=1;  j=2;
% first margin stuff
[theta1m,~,~,resids1] = ARMAX_2(rets1(:,i),mean_order(i,1),mean_order(i,2));
resids1 = [zeros(max(mean_order(i,1),mean_order(i,2)),1);resids1];
scores1m = LLgrad_1('ARMA_LLa',theta1m,mean_order(i,1),mean_order(i,2),rets1(:,i));
[theta1v, ~, hessinv1, robustSE1, hhat1, scores1v] = multigarch_AJP(resids1,1,1,1,'GJRGARCH','NORMAL');
hess1v = inv(hessinv1)/T;  % Kevin's code returns the inverse hessian - I convert it back to the hessian here, then divide by T
stdresids1 = resids1./sqrt(hhat1);    
theta1s = outSKEWT(i,:)';
theta1 = [theta1m;theta1v;theta1s];
scores1s = LLgrad_1('skewtdis_LLa',theta1s,stdresids1);
scores1 = [scores1m,scores1v,scores1s];
B1 = newey_west(scores1,floor(T^(1/3)));
H1 = zeros(size(scores1,2),size(scores1,2));
H1m = hessian('ARMA_LL',theta1m,mean_order(i,1),mean_order(i,2),rets1(:,i))/T;
H1v = hess1v;
H1s = hessian_2stage('skewtdis_ARMA_GJRGARCH_LL',theta1,length(theta1s),[],rets1(:,i),mean_order(i,1),mean_order(i,2),1,1,1)/T;
H1(1:length(theta1m),1:length(theta1m)) = H1m;
H1(length(theta1m)+1:length(theta1m)+length(theta1v)  ,length(theta1m)+1:length(theta1m)+length(theta1v)) = H1v;
H1(end-1:end,:) = H1s;
V1 = inv(H1)*B1*(inv(H1)')/T;  % avar[thetahat] = Ainv*B*Ainv, so V[thetahat] ~~ Ainv*B*Ainv/T
[theta1,sqrt(diag(V1)),theta1./sqrt(diag(V1))]

% second margin stuff
[theta2m,~,~,resids2] = ARMAX_2(rets1(:,j),mean_order(j,1),mean_order(j,2));
resids2 = [zeros(max(mean_order(j,1),mean_order(j,2)),1);resids2];
scores2m = LLgrad_1('ARMA_LLa',theta2m,mean_order(j,1),mean_order(j,2),rets1(:,j));
[theta2v, ~, hessinv2, robustSE2, hhat2, scores2v] = multigarch_AJP(resids2,1,1,1,'GJRGARCH','NORMAL');
hess2v = inv(hessinv2)/T;  % Kevin's code returns the inverse hessian - I convert it back to the hessian here, then divide by T
stdresids2 = resids2./sqrt(hhat2);    
theta2s = outSKEWT(j,:)';
theta2 = [theta2m;theta2v;theta2s];
scores2s = LLgrad_1('skewtdis_LLa',theta2s,stdresids2);
scores2 = [scores2m,scores2v,scores2s];
B2 = newey_west(scores2,floor(T^(1/3)));
H2 = zeros(size(scores2,2),size(scores2,2));
H2m = hessian('ARMA_LL',theta2m,mean_order(j,1),mean_order(j,2),rets1(:,j))/T;
H2v = hess2v;
H2s = hessian_2stage('skewtdis_ARMA_GJRGARCH_LL',theta2,length(theta2s),[],rets1(:,j),mean_order(j,1),mean_order(j,2),1,1,1)/T;
H2(1:length(theta2m),1:length(theta2m)) = H2m;
H2(length(theta2m)+1:length(theta2m)+length(theta2v)  ,length(theta2m)+1:length(theta2m)+length(theta2v)) = H2v;
H2(end-1:end,:) = H2s;
V2 = inv(H2)*B2*(inv(H2)')/T;  % avar[thetahat] = Ainv*B*Ainv, so V[thetahat] ~~ Ainv*B*Ainv/T
[theta2,sqrt(diag(V2)),theta2./sqrt(diag(V2))]


% copula stuff:
uu=2;  % parametric margins 

% Rot Gumbel copula
scoresc = LLgrad_1('Rotgumbel_GAS_CLa',kappa8b3,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
scoresALL = [scores1,scores2,scoresc];
BB = newey_west(scoresALL,floor(T^(1/3)));
thetac = kappa8b3;
kappa8b3'
tic;
Hc = hessian_2stage('margin_margin_copula_CL1',...
   [theta1;theta2;thetac],size(thetac,1),[],rets1(:,[i,j]),...
   'skewtdis_ARMA_GJRGARCH_LL','skewtdis_ARMA_GJRGARCH_LL','Rotgumbel_GAS_CL',...
   size(theta1,1),size(theta2,1),5,5,...
   mean_order(i,1),mean_order(i,2),1,1,1,...
   mean_order(j,1),mean_order(j,2),1,1,1,thetaALL(7,1,uu),[KK,HESSgumbel(:,1)])/T;
toc  % 53 seconds
thetaALLmle = [theta1;theta2;thetac];
HH = zeros(length(thetaALLmle),length(thetaALLmle));
HH(1:length(theta1),1:length(theta1)) = H1;
HH(length(theta1) + (1:length(theta2)),length(theta1) + (1:length(theta2))) = H2;
HH(end-length(thetac)+1:end,:) = Hc;
VALL = inv(HH)*BB*(inv(HH)')/T;
eig(VALL)
[theta1,sqrt(diag(V1)),sqrt(diag(VALL(1:length(theta1),1:length(theta1))))]
[theta2,sqrt(diag(V2)),sqrt(diag(VALL(1+length(theta1):length(theta1)+length(theta2),1+length(theta1):length(theta1)+length(theta2))))]
[thetaALLmle,sqrt(diag(VALL))]
table5([2;4;6],2) = sqrt(diag(VALL(end-2:end,end-2:end)));


% Student's t copula
% WARNING: THIS PART IS PRETTY SLOW
tic;scoresc = LLgrad_1('tcopula_GAS_CLa',theta8tv2,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);toc  % takes 35 seconds
scoresALL = [scores1,scores2,scoresc];
BB = newey_west(scoresALL,floor(T^(1/3)));
thetac = theta8tv2;  %  0.019872     0.065285      0.99121     0.088728
tic;
hstep0 = -eps.^(1/3)*max(abs([theta1;theta2;thetac]),1e-2);  
hstep5 = hstep0;
hstep5(end-length(thetac)+1:end) = -ones(length(thetac),1)*1e-3;  % it turned out to require some tinkering with the step size to get a reasonable-looking VCV matrix. This step size for the bottom rows works OK
Hc = hessian_2stage('margin_margin_copula_CL1',...
   [theta1;theta2;thetac],size(thetac,1),hstep5,rets1(:,[i,j]),...
   'skewtdis_ARMA_GJRGARCH_LL','skewtdis_ARMA_GJRGARCH_LL','tcopula_GAS_CL',...
   size(theta1,1),size(theta2,1),5,5,...
   mean_order(i,1),mean_order(i,2),1,1,1,...
   mean_order(j,1),mean_order(j,2),1,1,1,thetaALL(8,1,uu),RR,NN,HESSstudt)/T;toc  % takes 11 minutes

thetaALLmle = [theta1;theta2;thetac];
HH = zeros(length(thetaALLmle),length(thetaALLmle));
HH(1:length(theta1),1:length(theta1)) = H1;
HH(length(theta1) + (1:length(theta2)),length(theta1) + (1:length(theta2))) = H2;
HH(end-length(thetac)+1:end,:) = Hc;
VALL = inv(HH)*BB*(inv(HH)')/T;
eig(VALL)
table5([9;11;13;15],2) = sqrt(diag(Hc(end-3:end,end-3:end)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap std errors for param model 

uu=2;  % skew t margins

% WARNING: THIS PART IS VERY SLOW!
options = optimset('Display','iter','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'MaxIter',100,'DiffMaxChange',Inf,'DiffMinChange',0 );  % 24apr12: decreasing MaxIter here - it only goes above about 20 when there is a problem, so no need to keep it at default of 400.
bootreps = 100;  % should be at least 100
block_length = 60;  % need a longer block length here, given the persistence implied by the GAS parameters
bootdates = stat_bootstrap_function_21(T,bootreps,block_length);
kappa88boot = nan(bootreps,3+4+2);  % 24apr12: added two cols here to include the exit flags, to keep a check on how often there is a problem.
tic;
for bb=1:bootreps
    rets1boot = rets1(bootdates(:,bb),:);
    [~,~,~,resids_1boot] = ARMAX_2(rets1boot(:,1),2,0);
    resids_1boot = [zeros(2,1);resids_1boot];
    [~,~,~,resids_2boot] = ARMAX_2(rets1boot(:,2),0,0);
    [~, ~, ~, ~, hhat1boot] = multigarch_AJP(resids_1boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    [~, ~, ~, ~, hhat2boot] = multigarch_AJP(resids_2boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    stdresids1 = resids_1boot./sqrt(hhat1boot);
    stdresids2 = resids_2boot./sqrt(hhat2boot);
    
    lower = [2.1, -0.99]; upper = [Inf, 0.99 ];
    warning off;
    skewt1boot = fmincon('skewtdis_LL',outSKEWT(1,:)',[],[],[],[],lower,upper,[],options,stdresids1);
    warning off;
    skewt2boot = fmincon('skewtdis_LL',outSKEWT(2,:)',[],[],[],[],lower,upper,[],options,stdresids2);
    Uskewt1 = skewtdis_cdf(stdresids1,skewt1boot(1),skewt1boot(2));
    Uskewt2 = skewtdis_cdf(stdresids2,skewt2boot(1),skewt2boot(2));
    
    % 7. Rotated Gumbel GAS copula
    lower = [-3 ; -3 ; 0 ];
    upper = [ 3 ; 3 ; 0.999 ];
    [kappa8b3b,~,exitflag1] = fmincon('Rotgumbel_GAS_CL',kappa8b3,[],[],[],[],lower,upper,[],options,[Uskewt1,Uskewt2],thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
    
    % 8. Student's t copula
    lower = [-3 ; -3 ; 0     ; 0.01];
    upper = [ 3 ;  3 ; 0.999 ; 0.45];
    [theta8tv2b,~,exitflag2] = fmincon('tcopula_GAS_CL',theta8tv2,[],[],[],[],lower,upper,[],options,[Uskewt1,Uskewt2],thetaALL(8,1,uu),RR,NN,HESSstudt);
    
    kappa88boot(bb,:) = [kappa8b3b',theta8tv2b',exitflag1,exitflag2];
    
    [bb,toc]  
end
toc  % took 34 hours for bootreps=100.  SLOOOOOOOOOW!
table5([2;4;6],3) = std(kappa88boot(:,1:3));
table5([9;11;13;15],3) = std(kappa88boot(:,4:7));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation std errors for param model 

% WARNING: THIS PART IS VERY SLOW!

uu=2;  % skew t margins

bootreps = 100;
kappa88sim = nan(bootreps,3+4);
GOFsim100 = nan(bootreps,2,2);  % will also get info for GOF test that we'll use later (rather than run another slow simulation separately)
uu=2;
tic;
for bb=1:bootreps
    % First simulate from the copula
    U1 = Rotgumbel_GAS_CL(kappa7tvs,T,thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
    U2 = tcopula_GAS_rnd(theta8tvs,T,thetaALL(8,1,uu),RR,NN,HESSstudt);
    UUU = U1;
    UUU(:,:,2) = U2;
    % then obtain the std resids
    EEE = nan(size(UUU));  
    for cc=1:2;
        for mm=1:2;
            EEE(:,mm,cc) = skewtdis_inv(UUU(:,mm,cc),outSKEWT(mm,1),outSKEWT(mm,2));
        end
    end
    % next obtain the marginal dynamics
    MMM = nan(size(UUU));  % cond mean
    HHH = nan(size(UUU));  % vol
    YYY = nan(size(UUU));  % simulated raw data
    for cc=1:2
        for mm=1:2;
            HHH(1:2,mm,cc) = GARCHparams12(1,mm)/(1-GARCHparams12(2,mm)-GARCHparams12(3,mm)/2-GARCHparams12(4,mm));   % starting off at unconditional vol
            MMM(1:2,mm,cc) = mean(rets1(:,mm));  % starting off at unconditional mean
            YYY(1:2,mm,cc) = MMM(1:2,mm,cc) + sqrt(HHH(1:2,mm,cc)).*EEE(1:2,mm,cc);
            for tt=3:T
                HHH(tt,mm,cc) = GARCHparams12(1,mm) + GARCHparams12(2,mm)*HHH(tt-1,mm,cc)*(EEE(tt-1,mm,cc)^2) ...
                    + GARCHparams12(3,mm)*HHH(tt-1,mm,cc)*(EEE(tt-1,mm,cc)^2)*(EEE(tt-1,mm,cc)<0) ...
                    + GARCHparams12(4,mm)*HHH(tt-1,mm,cc);
                if mm==1  % then use AR(2) for mean
                    MMM(tt,mm,cc) = ARparams1(1) + ARparams1(2)*YYY(tt-1,mm,cc) + ARparams1(3)*YYY(tt-2,mm,cc);
                else
                    MMM(tt,mm,cc) = ARparams2(1);
                end
                YYY(tt,mm,cc) = MMM(tt,mm,cc) + sqrt(HHH(tt,mm,cc))*EEE(tt,mm,cc);
            end
        end
    end
    
    for cc=1:2
        % now estimate the models on the simulated data
        [~,~,~,resids_1boot] = ARMAX_2(YYY(:,1,cc),2,0);
        resids_1boot = [zeros(2,1);resids_1boot];
        [~,~,~,resids_2boot] = ARMAX_2(YYY(:,2,cc),0,0);
        [~, ~, ~, ~, hhat1boot] = multigarch_AJP(resids_1boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
        [~, ~, ~, ~, hhat2boot] = multigarch_AJP(resids_2boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
        stdresids1 = resids_1boot./sqrt(hhat1boot);
        stdresids2 = resids_2boot./sqrt(hhat2boot);
        
        lower = [2.1, -0.99]; upper = [Inf, 0.99 ];
        warning off;
        skewt1boot = fmincon('skewtdis_LL',outSKEWT(1,:)',[],[],[],[],lower,upper,[],options,stdresids1);
        warning off;
        skewt2boot = fmincon('skewtdis_LL',outSKEWT(2,:)',[],[],[],[],lower,upper,[],options,stdresids2);
        Uskewt1 = skewtdis_cdf(stdresids1,skewt1boot(1),skewt1boot(2));
        Uskewt2 = skewtdis_cdf(stdresids2,skewt2boot(1),skewt2boot(2));
    
        if cc==1
            % 7. Rotated Gumbel GAS copula
            lower = [-3 ; -3 ; 0 ];
            upper = [ 3 ; 3 ; 0.999 ];
            kappa8b3b = fmincon('Rotgumbel_GAS_CL',kappa7tvs,[],[],[],[],lower,upper,[],options,[Uskewt1,Uskewt2],thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
            
            % getting the Rosenblatt transformations to use in GOF tests
            [CL,kappat,ft] = Rotgumbel_GAS_CL(kappa8b3b,[Uskewt1,Uskewt2],thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);  % getting the time path of the copula parameter
            Vhat1 = 1-Uskewt1;
            Vhat2 = GumbelUgivenV_t(1-Uskewt2,1-Uskewt1,0,kappat);  % conditional distribution of U2 | U1
        elseif cc==2
            % 8. Student's t copula
            lower = [-3 ; -3 ; 0     ; 0.01];
            upper = [ 3 ;  3 ; 0.999 ; 0.45];
            theta8tv2b = fmincon('tcopula_GAS_CL',theta8tvs,[],[],[],[],lower,upper,[],options,[Uskewt1,Uskewt2],thetaALL(8,1,uu),RR,NN,HESSstudt);
            
            % getting the Rosenblatt transformations to use in GOF tests
            [CL,rhot,ft] = tcopula_GAS_CL(theta8tv2b,[Uskewt1,Uskewt2],thetaALL(8,1,uu),RR,NN,HESSstudt);
            Vhat1 = Uskewt1;
            Vhat2 = nan(T,1);
            for tt=1:T;  % need to loop over t as my "mvt_cond_cdf" functino does not accept time-varying inputs (it assumes that the correlatino and DoF are constant)s
                Vhat2(tt) = mvt_cond_cdf(tinv(Uskewt2(tt),1/theta8tv2b(end)),tinv(Uskewt1(tt),1/theta8tv2b(end)),zeros(1,2),[[1,rhot(tt)];[rhot(tt),1]],1/theta8tv2b(end));
            end
        end
        [KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2],'IndepCop_cdf');
        GOFsim100(bb,:,cc) = [KSstat1 CVMstat1];
    end
    kappa88sim(bb,:) = [kappa8b3b',theta8tv2b'];
    
    [bb,toc]    % about ?? mins per replication
    
end
toc  % took 27.8 hours for bootreps=100

table5([2;4;6],4) = std(kappa88sim(:,1:3));
table5([9;11;13;15],4) = std(kappa88sim(:,4:7));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive std errors for semiparam model 

uu=1;  % EDF margins

% 7. Rot Gumbel copula
scoresc = LLgrad_1('Rotgumbel_GAS_CLa',kappa7tvs,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
Hc1 =     hessian('Rotgumbel_GAS_CL',  kappa7tvs,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)])/T;  % used in naive VCV matrix for copula params
Bc = newey_west(scoresc,floor(T^(1/3)));
Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
[kappa7tvs, sqrt(diag(Vc2)) ]
table5([2;4;6],5) = sqrt(diag(Vc2));


% 8. Student's t copula
scoresc = LLgrad_1('tcopula_GAS_CLa',theta8tvs,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);
Hc1 =      hessian('tcopula_GAS_CL', theta8tvs,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt)/T;  % used in naive VCV matrix for copula params. NOTE: STEP SIZE OF h = -0.0001; WAS THE ONE THAT WORKED HERE (AFTER ABOUT 4 TRIES)
Bc = newey_west(scoresc,floor(T^(1/3)));
Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
[theta8tvs, sqrt(diag(Vc2))]
table5([9;11;13;15],5) = sqrt(diag(Vc2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bootstrap std errors for semiparam model 

% WARNING: THIS PART IS VERY SLOW!

uu=1;  % EDF margins

bootreps = 100;  % should be at least 100
block_length = 60;
bootdates = stat_bootstrap_function_21(T,bootreps,block_length);
kappa88Sboot = nan(bootreps,3+4);
tic;
for bb=1:bootreps
    rets1boot = rets1(bootdates(:,bb),:);
    [~,~,~,resids_1boot] = ARMAX_2(rets1boot(:,1),2,0);
    resids_1boot = [zeros(2,1);resids_1boot];
    [~,~,~,resids_2boot] = ARMAX_2(rets1boot(:,2),0,0);
    [~, ~, ~, ~, hhat1boot] = multigarch_AJP(resids_1boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    [~, ~, ~, ~, hhat2boot] = multigarch_AJP(resids_2boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    stdresids1 = resids_1boot./sqrt(hhat1boot);
    stdresids2 = resids_2boot./sqrt(hhat2boot);
    
    UbootS = empiricalCDF([stdresids1,stdresids2]);  
    
    % 7. Rotated Gumbel GAS copula
    lower = [-3 ; -3 ; 0 ];
    upper = [ 3 ; 3 ; 0.999 ];
    kappa8b3b = fmincon('Rotgumbel_GAS_CL',kappa7tvs,[],[],[],[],lower,upper,[],options,UbootS,thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
    
    % 8. Student's t copula
    lower = [-3 ; -3 ; 0     ; 0.01];
    upper = [ 3 ;  3 ; 0.999 ; 0.45];
    theta8tv2b = fmincon('tcopula_GAS_CL',theta8tvs,[],[],[],[],lower,upper,[],options,UbootS,thetaALL(8,1,uu),RR,NN,HESSstudt);
    
    kappa88Sboot(bb,:) = [kappa8b3b',theta8tv2b'];
    
    [bb,toc]    % about 15.7 mins per replication
end
toc  % took 27.8 hours for 100 reps. SLOW!!

table5([2;4;6],6) = std(kappa88Sboot(:,1:3));
table5([9;11;13;15],6) = std(kappa88Sboot(:,4:7));


% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_9.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_9.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GOODNESS OF FIT TESTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table6 = nan(6+6,4);  % first 6 rows are parametric, last 6 are semipara ; cols are KSc, CVMc, KSr, CVMr
% Note that I'm not replicating the "naive" test results reported in the chapter here, as those also require a simulation but are not correct 
% (so they're a pain AND they're wrong!)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametric margins, constant copula
% KS and CVM on the copula directly

uu=2;  % skew t margins
[KSstat1 CVMstat1] = copula_GOF_stats(Uall(:,:,uu),'NormalCopula_cdf',thetaALL(1,1,uu));
[KSstat2 CVMstat2] = copula_GOF_stats(Uall(:,:,uu),'clayton_cdf',thetaALL(2,1,uu));
[KSstat3 CVMstat3] = copula_GOF_stats(1-Uall(:,:,uu),'gumbel_cdf',thetaALL(7,1,uu));
[KSstat4 CVMstat4] = copula_GOF_stats(Uall(:,:,uu),'tCopula_cdf_new',thetaALL(8,1:2,uu)');

GOFstats111 = [[KSstat1;KSstat2;KSstat3;KSstat4],[CVMstat1;CVMstat2;CVMstat3;CVMstat4]];
for mm=1:4;
    for tt=1:2;
        table6(mm,tt) = mean(GOFsim1(:,tt,mm)>GOFstats111(mm,tt));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonparametric margins, constant copula
% KS and CVM on the copula directly

uu=1; % EDF margins
[KSstat1 CVMstat1] = copula_GOF_stats(Uall(:,:,uu),'NormalCopula_cdf',thetaALL(1,1,uu));
[KSstat2 CVMstat2] = copula_GOF_stats(Uall(:,:,uu),'clayton_cdf',thetaALL(2,1,uu));
[KSstat3 CVMstat3] = copula_GOF_stats(1-Uall(:,:,uu),'gumbel_cdf',thetaALL(7,1,uu));
[KSstat4 CVMstat4] = copula_GOF_stats(Uall(:,:,uu),'tCopula_cdf_new',thetaALL(8,1:2,uu)');

GOFstats111 = [[KSstat1;KSstat2;KSstat3;KSstat4],[CVMstat1;CVMstat2;CVMstat3;CVMstat4]];
for mm=1:4;
    for tt=1:2;
        table6(6+mm,tt) = mean(GOFsim3(:,tt,mm)>GOFstats111(mm,tt));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametric margins, constant copula
% KS and CVM on the Rosenblatt transforms

uu=2; % skew t margins
Vhat1 = Uall(:,1,uu);
Vhat2all = nan(T,4);
for cc=1:4;
    if cc==1;
        Vhat2all(:,cc) = normcdf(norminv(Uall(:,2,uu)),thetaALL(1,1,uu)*norminv(Uall(:,1,uu)),sqrt(1-thetaALL(1,1,uu)^2));  % conditional copula of V2 | V1
    elseif cc==2
        Vhat2all(:,cc) = ClaytonUgivenV_t(Uall(:,2,uu),Uall(:,1,uu),0,thetaALL(2,1,uu));  % conditional distribution of U2 | U1
    elseif cc==3
        Vhat2all(:,cc) = GumbelUgivenV_t(1-Uall(:,2,uu),1-Uall(:,1,uu),0,thetaALL(7,1,uu));  % conditional distribution of U2 | U1
    elseif cc==4;
        Vhat2all(:,cc) = mvt_cond_cdf(tinv(Uall(:,2,uu),1/thetaALL(8,2,uu)),tinv(Uall(:,1,uu),1/thetaALL(8,2,uu)),zeros(1,2),[[1,thetaALL(8,1,uu)];[thetaALL(8,1,uu),1]],1/thetaALL(8,2,uu));
    end
end

[KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2all(:,1)],'IndepCop_cdf');
[KSstat2 CVMstat2] = copula_GOF_stats([Vhat1,Vhat2all(:,2)],'IndepCop_cdf');
[KSstat3 CVMstat3] = copula_GOF_stats([Vhat1,Vhat2all(:,3)],'IndepCop_cdf');
[KSstat4 CVMstat4] = copula_GOF_stats([Vhat1,Vhat2all(:,4)],'IndepCop_cdf');

GOFstats111 = [[KSstat1;KSstat2;KSstat3;KSstat4],[CVMstat1;CVMstat2;CVMstat3;CVMstat4]];
for mm=1:4;
    for tt=1:2;
        table6(mm,2+tt) = mean(GOFsim10(:,tt,mm)>GOFstats111(mm,tt));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonparametric margins, constant copula
% KS and CVM on the Rosenblatt transforms

uu=1; % EDF margins
Vhat1 = Uall(:,1,uu);
Vhat2all = nan(T,4);
for cc=1:4;
    if cc==1;
        Vhat2all(:,cc) = normcdf(norminv(Uall(:,2,uu)),thetaALL(1,1,uu)*norminv(Uall(:,1,uu)),sqrt(1-thetaALL(1,1,uu)^2));  % conditional copula of V2 | V1
    elseif cc==2
        Vhat2all(:,cc) = ClaytonUgivenV_t(Uall(:,2,uu),Uall(:,1,uu),0,thetaALL(2,1,uu));  % conditional distribution of U2 | U1
    elseif cc==3
        Vhat2all(:,cc) = GumbelUgivenV_t(1-Uall(:,2,uu),1-Uall(:,1,uu),0,thetaALL(7,1,uu));  % conditional distribution of U2 | U1
    elseif cc==4;
        Vhat2all(:,cc) = mvt_cond_cdf(tinv(Uall(:,2,uu),1/thetaALL(8,2,uu)),tinv(Uall(:,1,uu),1/thetaALL(8,2,uu)),zeros(1,2),[[1,thetaALL(8,1,uu)];[thetaALL(8,1,uu),1]],1/thetaALL(8,2,uu));
    end
end

[KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2all(:,1)],'IndepCop_cdf');
[KSstat2 CVMstat2] = copula_GOF_stats([Vhat1,Vhat2all(:,2)],'IndepCop_cdf');
[KSstat3 CVMstat3] = copula_GOF_stats([Vhat1,Vhat2all(:,3)],'IndepCop_cdf');
[KSstat4 CVMstat4] = copula_GOF_stats([Vhat1,Vhat2all(:,4)],'IndepCop_cdf');

GOFstats111 = [[KSstat1;KSstat2;KSstat3;KSstat4],[CVMstat1;CVMstat2;CVMstat3;CVMstat4]];
for mm=1:4;
    for tt=1:2;
        table6(6+mm,2+tt) = mean(GOFsim30(:,tt,mm)>GOFstats111(mm,tt));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametric margins, time varying copula
% KS and CVM on the Rosenblatt transforms

uu=2;  % skew t margins
Vhat1 = Uall(:,1,uu);
Vhat2all = nan(T,4);
for cc=1:2;
    if cc==1
        [CL,kappat,ft] = Rotgumbel_GAS_CL(kappa8b3,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);  % getting the time path of the copula parameter
        Vhat2all(:,cc) = GumbelUgivenV_t(1-Uall(:,2,cc),1-Uall(:,1,cc),0,kappat);  % conditional distribution of U2 | U1
    elseif cc==2;
        [CL,rhot,ft] = tcopula_GAS_CL(theta8tv2,Uall(:,:,cc),thetaALL(8,1,uu),RR,NN,HESSstudt);
        for tt=1:T;  % need to loop over t as my "mvt_cond_cdf" functino does not accept time-varying inputs (it assumes that the correlatino and DoF are constant)s
            Vhat2all(tt,cc) = mvt_cond_cdf(tinv(Uall(tt,2,uu),1/theta8tv2(end)),tinv(Uall(tt,1,uu),1/theta8tv2(end)),zeros(1,2),[[1,rhot(tt)];[rhot(tt),1]],1/theta8tv2(end));
        end
    end
end

[KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2all(:,1)],'IndepCop_cdf');
[KSstat2 CVMstat2] = copula_GOF_stats([Vhat1,Vhat2all(:,2)],'IndepCop_cdf');

GOFstats111 = [[KSstat1;KSstat2],[CVMstat1;CVMstat2]];
for mm=1:2;
    for tt=1:2;
        table6(4+mm,2+tt) = mean(GOFsim100(:,tt,mm)>GOFstats111(mm,tt));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonparametric margins, time varying copula
% KS and CVM on the Rosenblatt transforms

% Note that the method below does not (yet) have theoretical foundation. I am treating the EDF below in the same way as the parametric margin case,
% but this approach still needs a formal result showing that it is appropriate.

uu=1; % EDF margins

% WARNING: THIS PART IS VERY SLOW!

bootreps = 100;
GOFsim120 = nan(bootreps,2,2);
tic;
for bb=1:bootreps
    % First simulate from the copula
    U1 = Rotgumbel_GAS_CL(kappa7tvs,T,thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
    U2 = tcopula_GAS_rnd(theta8tvs,T,thetaALL(8,1,uu),RR,NN,HESSstudt);
    UUU = U1;
    UUU(:,:,2) = U2;
    % then obtain the std resids
    EEE = nan(size(UUU));  
    for cc=1:2;
        for mm=1:2;
            EEE(:,mm,cc) = quantile(stdresids(:,mm),UUU(:,mm,cc),1);  % using empirical dist fn to get std resids (and using linear interpolation to smooth it out a bit). This is where the nonparam bit is
        end
    end
    % next obtain the marginal dynamics
    MMM = nan(size(UUU));  % cond mean
    HHH = nan(size(UUU));  % vol
    YYY = nan(size(UUU));  % simulated raw data
    for cc=1:2
        for mm=1:2;
            HHH(1:2,mm,cc) = GARCHparams12(1,mm)/(1-GARCHparams12(2,mm)-GARCHparams12(3,mm)/2-GARCHparams12(4,mm));   % starting off at unconditional vol
            MMM(1:2,mm,cc) = mean(rets1(:,mm));  % starting off at unconditional mean
            YYY(1:2,mm,cc) = MMM(1:2,mm,cc) + sqrt(HHH(1:2,mm,cc)).*EEE(1:2,mm,cc);
            for tt=3:T
                HHH(tt,mm,cc) = GARCHparams12(1,mm) + GARCHparams12(2,mm)*HHH(tt-1,mm,cc)*(EEE(tt-1,mm,cc)^2) ...
                    + GARCHparams12(3,mm)*HHH(tt-1,mm,cc)*(EEE(tt-1,mm,cc)^2)*(EEE(tt-1,mm,cc)<0) ...
                    + GARCHparams12(4,mm)*HHH(tt-1,mm,cc);
                if mm==1  % then use AR(2) for mean
                    MMM(tt,mm,cc) = ARparams1(1) + ARparams1(2)*YYY(tt-1,mm,cc) + ARparams1(3)*YYY(tt-2,mm,cc);
                else
                    MMM(tt,mm,cc) = ARparams2(1);
                end
                YYY(tt,mm,cc) = MMM(tt,mm,cc) + sqrt(HHH(tt,mm,cc))*EEE(tt,mm,cc);
            end
        end
    end
    
    % now estimate the models on the simulated data
    [~,~,~,resids_1boot] = ARMAX_2(YYY(:,1,cc),2,0);
    resids_1boot = [zeros(2,1);resids_1boot];
    [~,~,~,resids_2boot] = ARMAX_2(YYY(:,2,cc),0,0);
    [~, ~, ~, ~, hhat1boot] = multigarch_AJP(resids_1boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    [~, ~, ~, ~, hhat2boot] = multigarch_AJP(resids_2boot,1,1,1,'GJRGARCH','NORMAL',[],parameters4);
    stdresids1 = resids_1boot./sqrt(hhat1boot);
    stdresids2 = resids_2boot./sqrt(hhat2boot);

    Uhat = empiricalCDF([stdresids1,stdresids2]);  % using EDF to get unif variables
    
    for cc=1:2
        if cc==1
            % 7. Rotated Gumbel GAS copula
            lower = [-3 ; -3 ; 0 ];
            upper = [ 3 ; 3 ; 0.999 ];
            kappa8b3b = fmincon('Rotgumbel_GAS_CL',kappa7tvs,[],[],[],[],lower,upper,[],options,Uhat,thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
            [CL,kappat,ft] = Rotgumbel_GAS_CL(kappa8b3b,Uhat,thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);  % getting the time path of the copula parameter

            % getting the Rosenblatt transformations
            Vhat1 = 1-Uhat(:,1);
            Vhat2 = GumbelUgivenV_t(1-Uhat(:,2),1-Uhat(:,1),0,kappat);  % conditional distribution of U2 | U1
        elseif cc==2
            % 8. Student's t copula
            lower = [-0 ; -0.01  ; 0     ; 0.01];
            upper = [ 0.25 ;  0.25 ; 0.999 ; 0.45];
            theta8tv2b = fmincon('tcopula_GAS_CL',theta8tvs,[],[],[],[],lower,upper,[],options,Uhat,thetaALL(8,1,uu),RR,NN,HESSstudt);
            [CL,rhot,ft] = tcopula_GAS_CL(theta8tv2b,Uhat,thetaALL(8,1,uu),RR,NN,HESSstudt);
            
            % getting the Rosenblatt transformations
            Vhat1 = Uhat(:,1);
            Vhat2 = nan(T,1);
            for tt=1:T;  % need to loop over t as my "mvt_cond_cdf" functino does not accept time-varying inputs (it assumes that the correlatino and DoF are constant)s
                Vhat2(tt) = mvt_cond_cdf(tinv(Uhat(tt,2),1/theta8tv2b(end)),tinv(Uhat(tt,1),1/theta8tv2b(end)),zeros(1,2),[[1,rhot(tt)];[rhot(tt),1]],1/theta8tv2b(end));
            end
        end
        [KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2],'IndepCop_cdf');
        GOFsim120(bb,:,cc) = [KSstat1 CVMstat1];
    end
    [bb,toc]    % about 29 mins per replication
end
toc  % took 27.2 hours

uu=1;

Vhat1 = Uall(:,1,uu);
Vhat2all = nan(T,4);
for cc=1:2;
    if cc==1
        [CL,kappat,ft] = Rotgumbel_GAS_CL(kappa7tvs,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);  % getting the time path of the copula parameter
        Vhat2all(:,cc) = GumbelUgivenV_t(1-Uall(:,2,cc),1-Uall(:,1,cc),0,kappat);  % conditional distribution of U2 | U1
    elseif cc==2;
        [CL,rhot,ft] = tcopula_GAS_CL(theta8tvs,Uall(:,:,cc),thetaALL(8,1,uu),RR,NN,HESSstudt);
        for tt=1:T;  % need to loop over t as my "mvt_cond_cdf" functino does not accept time-varying inputs (it assumes that the correlatino and DoF are constant)s
            Vhat2all(tt,cc) = mvt_cond_cdf(tinv(Uall(tt,2,uu),1/theta8tv2(end)),tinv(Uall(tt,1,uu),1/theta8tv2(end)),zeros(1,2),[[1,rhot(tt)];[rhot(tt),1]],1/theta8tv2(end));
        end
    end
end

[KSstat1 CVMstat1] = copula_GOF_stats([Vhat1,Vhat2all(:,1)],'IndepCop_cdf');
[KSstat2 CVMstat2] = copula_GOF_stats([Vhat1,Vhat2all(:,2)],'IndepCop_cdf');

GOFstats111 = [[KSstat1;KSstat2],[CVMstat1;CVMstat2]];
for mm=1:2;
    for tt=1:2;
        table6(6+4+mm,2+tt) = mean(GOFsim120(:,tt,mm)>GOFstats111(mm,tt));
    end
end

info.fmt = strvcat('%10.2f');
info.cnames = strvcat('KS_C','CVM_C','KS_R','CVM_R');
info.rnames = strvcat('.','Normal','Clayton','Rot Gumbel','Stud t','RotGum-GAS','Studt-GAS');
fprintf(['\n\nTable 6a: Goodness of fit tests for copula models - parametric margins\n'])
mprint(table6(1:6,:),info)
fprintf(['\n\nTable 6b: Goodness of fit tests for copula models - nonparametric margins\n'])
mprint(table6(7:12,:),info)


% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_10.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_10.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IN SAMPLE MODEL COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table7 = nan(6+2+6+2,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametric margins

uu=2;  % parametric margins: Rivers and Vuong test (simple t test)

LL1a = NormalCopula_CLa(thetaALL(1,1,uu),Uall(:,:,uu));
LL2a = claytonCLa(thetaALL(2,1,uu),Uall(:,:,uu));
LL7a = gumbelCLa(thetaALL(7,1,uu),1-Uall(:,:,uu));
LL8a = tcopulaCL2a(thetaALL(8,1:2,uu),Uall(:,:,uu));

LL70a = Rotgumbel_GAS_CLa(kappa8b3,Uall(:,:,uu),thetaALL(7,1,uu),[KK,HESSgumbel(:,1)]);
LL80a = tcopula_GAS_CLa(theta8tv2,Uall(:,:,uu),thetaALL(8,1,uu),RR,NN,HESSstudt);

LLALLa = -[LL1a,LL2a,LL7a,LL8a,LL70a,LL80a];

table7a = nan(6,6);
for jj=1:5;
    for ii=jj+1:6;
        temp = nwest(LLALLa(:,ii)-LLALLa(:,jj),ones(length(LLALLa),1),10);
        table7a(ii,jj) = temp.tstat;
    end
end

table7(1:4,1:4) = table7a(1:4,1:4);
table7(5,1:4) = sum(LLALLa(:,1:4));
table7(6,1:4) = ranks(-sum(LLALLa(:,1:4))');
table7(4,1) = thetaALL(8,2,uu)/table4(13,2);  % comparisno of Student's t with Normal is done via a simple t test on the "nu_inverse" parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonparametric margins
% Chen and Fan test (need more complicated estimate of vairance of tstat)

uu=1;  % EDF margins

LL1a = NormalCopula_CLa(thetaALL(1,1,uu),Uall(:,:,uu));
LL2a = claytonCLa(thetaALL(2,1,uu),Uall(:,:,uu));
LL7a = gumbelCLa(thetaALL(7,1,uu),1-Uall(:,:,uu));
LL8a = tcopulaCL2a(thetaALL(8,1:2,uu),Uall(:,:,uu));

LLALLa = -[LL1a,LL2a,LL7a,LL8a];
table7(6+5,1:4) = sum(LLALLa(:,1:4));
table7(6+6,1:4) = ranks(-sum(LLALLa(:,1:4))');

temp124 = [1;2;7;8];
tic;
ii=1;jj=2;[teststats,pvals] = chen_fan_PLR_test(thetaALL(temp124(ii),1,uu),thetaALL(temp124(jj),1,uu),Uall(:,:,uu),'NormalCopula_CLa','claytonCLa');table7(6+jj,ii) = teststats(2);
ii=1;jj=3;[teststats,pvals] = chen_fan_PLR_test(thetaALL(temp124(ii),1,uu),thetaALL(temp124(jj),1,uu),Uall(:,:,uu),'NormalCopula_CLa','Rotgumbel_GAS_CLa');table7(6+jj,ii) = teststats(2);
ii=2;jj=3;[teststats,pvals] = chen_fan_PLR_test(thetaALL(temp124(ii),1,uu),thetaALL(temp124(jj),1,uu),Uall(:,:,uu),'claytonCLa','Rotgumbel_GAS_CLa');table7(6+jj,ii) = teststats(2);

ii=1;jj=4;[teststats,pvals] = chen_fan_PLR_test(thetaALL(temp124(ii),1,uu),thetaALL(temp124(jj),1:2,uu)',Uall(:,:,uu),'NormalCopula_CLa','tcopulaCL2a');table7(6+jj,ii) = teststats(2);
ii=2;jj=4;[teststats,pvals] = chen_fan_PLR_test(thetaALL(temp124(ii),1,uu),thetaALL(temp124(jj),1:2,uu)',Uall(:,:,uu),'claytonCLa','tcopulaCL2a');table7(6+jj,ii) = teststats(2);
ii=3;jj=4;[teststats,pvals] = chen_fan_PLR_test(thetaALL(temp124(ii),1,uu),thetaALL(temp124(jj),1:2,uu)',Uall(:,:,uu),'Rotgumbel_GAS_CLa','tcopulaCL2a');table7(6+jj,ii) = teststats(2);
toc  % all six tests take about 45 seconds

info.fmt = strvcat('%10.2f');
info.cnames = strvcat('Normal','Clayton','Rot Gumbel','Stud t');
info.rnames = strvcat('.','Normal','Clayton','Rot Gumbel','Stud t','log L','Rank');
fprintf(['\n\nTable 7a: In-sample model comparisons for constant copula models - parametric margins\n'])
mprint(table7(1:6,:),info)
fprintf(['\n\nTable 7b: In-sample model comparisons for constant copula models - nonparametric margins\n'])
mprint(table7(7:12,:),info)

fprintf(['\n\nTable 7c: In-sample model comparison of Rotated Gumbel and Student''s t GAS copula models - parametric margins\n'])
fprintf(['t-statistic = ',num2str(table7a(6,5)),'\n'])


% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_11.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_11.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT OF SAMPLE MODEL COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table8 = nan(6+2+6+2+1,6);

min(dates1)  % 19950817
max(dates1)  % 20110520

% let's use first 10 years for estimation and remainder (about 5 years) for evaluation
Rsample = (dates1<=20050817);
Psample = (dates1>20050817);
R = sum(Rsample)  % 2519
P = sum(Psample)  % 1450
% looks like a good split

% I will use the same models as I considered for the full sample, though in practice one should really go back and find the optimal mean and variance
% models using just the in-sample period

% For simplicity I will consider a "fixed window" estimation strategy, so I will estimate the models just *once*, and then use those parameter
% estimates throughout the OOS period

% mean models
temp1 = ols(rets1(3:R,1),[ones(R-2,1),rets1(2:R-1,1),rets1(1:R-2,1)]);
MEANparams1 = temp1.beta;
MEANparams2 =mean(rets1(1:R,2));
residsIS = [[0;0;temp1.resid],rets1(1:R,2)-MEANparams2];

% vol models
ii=1;[GARCHparams1IS, ~, ~, ~, hhat1IS] = multigarch_AJP(residsIS(:,ii),1,1,1,'GJRGARCH','NORMAL',[],GARCHparams1);
ii=2;[GARCHparams2IS, ~, ~, ~, hhat2IS] = multigarch_AJP(residsIS(:,ii),1,1,1,'GJRGARCH','NORMAL',[],GARCHparams2);
stdresidsIS = residsIS./sqrt([hhat1IS,hhat2IS]);
GARCHparamsIS = [GARCHparams1IS,GARCHparams2IS];

% skewt marginal models
lower = [2.1, -0.99];
upper = [Inf, 0.99 ];
ii=1;theta1IS = fmincon('skewtdis_LL',outSKEWT(ii,:)',[],[],[],[],lower,upper,[],options,stdresidsIS(:,ii));
ii=2;theta2IS = fmincon('skewtdis_LL',outSKEWT(ii,:)',[],[],[],[],lower,upper,[],options,stdresidsIS(:,ii));
outSKEWT
outSKEWTIS = [theta1IS ,theta2IS ]'

% PIT variables
UhatIS = nan(T,2,2);  % day ; [margin1, margin2] ; [nonparam , skewt]
for mm=1:2;
    UhatIS(1:R,mm,1) = empiricalCDF(stdresidsIS(:,mm));
    UhatIS(1:R,mm,2) = skewtdis_cdf(stdresidsIS(:,mm),outSKEWTIS(mm,1),outSKEWTIS(mm,2));
end

% copula models
KAPPAhatIS = nan(4,6,2);  % parameter (nan if num params<4) ; [Normal, Clayton, RGum, Studt, RGum-GAS, Studt-GAS] ; [nonparam, skewt]
tic;
for uu=1:2
    % 1. Normal Copula
    KAPPAhatIS(1,1,uu) = corrcoef12(norminv(UhatIS(1:R,:,uu)));
    % 2. Clayton's copula
    lower = 0.0001;    warning off;
    KAPPAhatIS(1,2,uu) = fmincon('claytonCL',thetaALL(2,1,2),[],[],[],[],lower,[],[],options,UhatIS(1:R,:,uu));
    % 7. Rotated Gumbel copula
    lower = 1.1;
    upper = 5;
    KAPPAhatIS(1,3,uu) = fmincon('gumbelCL',thetaALL(7,1,2),[],[],[],[],lower,upper,[],options,1-UhatIS(1:R,:,uu));
    % 8. Student's t copula
    lower = [-0.9 , 0.01 ];
    upper = [ 0.9 , 0.45 ];
    KAPPAhatIS(1:2,4,uu) = fmincon('tcopulaCL2',thetaALL(8,:,2)',[],[],[],[],lower,upper,[],options,UhatIS(1:R,:,uu));
    
    % 7. Rotated Gumbel GAS copula
    lower = [-3 ; -3 ; 0 ];
    upper = [ 3 ; 3 ; 0.999 ];
    KAPPAhatIS(1:3,5,uu) = fmincon('Rotgumbel_GAS_CL',kappa8b3,[],[],[],[],lower,upper,[],options,UhatIS(1:R,:,uu),KAPPAhatIS(1,3,uu),[KK,HESSgumbel(:,1)]);
    % 8. Student's t copula
    lower = [-3 ; -3 ; 0     ; 0.01];
    upper = [ 3 ;  3 ; 0.999 ; 0.45];
    KAPPAhatIS(1:4,6,uu) = fmincon('tcopula_GAS_CL',theta8tv2,[],[],[],[],lower,upper,[],options,UhatIS(1:R,:,uu),KAPPAhatIS(1,4,uu),RR,NN,HESSstudt);
end
toc  % takes 17 mins for all 6 models and both types of margins (param, nonparm)
KAPPAhatIS

% OK, now we use these estimated parameters to obtain estimated mean, variance, and copula parameters through the OOS period
MUhatOOS = nan(T,2);
VOL2hatOOS = nan(T,2);
KAPPAhatOOS = nan(T,2,6,2);  % day ; [param1 (eg, kappa or rho), param2 (nu, for Stud t)] ; [Normal, Clayton, RGum, Studt, RGum-GAS, Studt-GAS] ;  [nonparam, skewt]

MUhatOOS(:,2) = MEANparams2;  % second margin had constant conditional mean
MUhatOOS(:,1) = [rets1(1:2,1); [ones(T-2,1),rets1(2:end-1,1),rets1(1:end-2,1)]*MEANparams1 ];  % cond mean for first margin
residsOOS = rets1 - MUhatOOS;

VOL2hatOOS(1:R,:) = [hhat1IS,hhat2IS];  
for mm=1:2;
    for tt=R+1:T;
        VOL2hatOOS(tt,mm) = GARCHparamsIS(1,mm) + GARCHparamsIS(2,mm)*(residsOOS(tt-1,mm)^2) ...
            + GARCHparamsIS(3,mm)*(residsOOS(tt-1,mm)^2)*(residsOOS(tt-1,mm)<0) ...
            + GARCHparamsIS(4,mm)*VOL2hatOOS(tt-1,mm);
    end
end
stdresidsOOS = residsOOS./sqrt(VOL2hatOOS);

UhatOOS = nan(T,2,2);
for mm=1:2;
    UhatOOS(1:R,mm,:) = UhatIS(1:R,mm,:);
    
    UhatOOS(R+1:end,mm,1) = empiricalCDF(stdresidsOOS(1:R,mm),stdresidsOOS(R+1:end,mm));  % EDF from in-sample data evaluated at OOS std resids
    UhatOOS(R+1:end,mm,2) = skewtdis_cdf(stdresidsOOS(R+1:end,mm),outSKEWTIS(mm,1),outSKEWTIS(mm,2));  % skew t dis with IS params evaluated at OOS std resids
end


% when using in-sample EDF, there is a chance that an OOS PIT observation gets a value of 0, which causes problems below. so check for this:
min(UhatOOS)  % so there may be such a problem with the EDF margins, but OK for param margins

uu=1;
mm=1;    
temp124 = find(UhatOOS(R+1:end,mm,uu)==0);
size(temp124)% just one obs
temp124  % 383
UhatOOS(R+temp124,mm,uu)
UhatOOS(R+temp124,mm,uu) = 1/(2*R);  % pushing this below 1/R, but above 0. This maintains its rank

mm=2;    
temp124 = find(UhatOOS(R+1:end,mm,uu)==0);
size(temp124)% just one obs
temp124  % again 383
UhatOOS(R+temp124,mm,uu)
UhatOOS(R+temp124,mm,uu) = 1/(2*R);  % pushing this below 1/R. This maintains its rank

% now getting parameters of the conditional copula
tic;
for uu=1:2;
    KAPPAhatOOS(:,1,1:4,uu) = ones(T,1)*KAPPAhatIS(1,1:4,uu);  % first four copulas all have constant parameters
    KAPPAhatOOS(:,2,4,uu)   = ones(T,1)*KAPPAhatIS(2,  4,uu);  % fourth copula (stud t) has constant DoF parameter
    KAPPAhatOOS(:,2,6,uu)   = ones(T,1)*KAPPAhatIS(4,  6,uu);  % sixth copula (stud t-GAS) has constant DoF parameter
    
    % now need to get the time-varying copula parameters for the RGum and Studt GAS models
    [~,kappatIS] = Rotgumbel_GAS_CL(KAPPAhatIS(1:3,5,uu),UhatOOS(1:R,:,uu),KAPPAhatIS(1,3,uu),[KK,HESSgumbel(:,1)]);
    KAPPAhatOOS(1:R,1,5,uu) = kappatIS;
    [~,rhotIS] = tcopula_GAS_CL( KAPPAhatIS(1:4,6,uu),UhatOOS(1:R,:,uu),KAPPAhatIS(1,4,uu),RR,NN,HESSstudt);
    KAPPAhatOOS(1:R,1,6,uu) = rhotIS;
    for tt=R+1:T;
        kappatOOS = rotgumbel_GAS_one_step(KAPPAhatIS(1:3,5,uu),KAPPAhatOOS(tt-1,1,5,uu),UhatOOS(tt-1,:,uu),KAPPAhatIS(1,3,uu),[KK,HESSgumbel(:,1)]);
        KAPPAhatOOS(tt,1,5,uu) = kappatOOS;
        rhotOOS = tcopula_GAS_one_step( KAPPAhatIS(1:4,6,uu),KAPPAhatOOS(tt-1,1,6,uu),UhatOOS(tt-1,:,uu),KAPPAhatIS(1,4,uu),RR,NN,HESSstudt);
        KAPPAhatOOS(tt,1,6,uu) = rhotOOS;
    end
end
toc  % takes 15 seconds 

figure(20201),plot([kappat,squeeze(KAPPAhatOOS(:,1,[3,5],1))]),legend('In sample GAS','OOS const','OOS GAS',2),title('Rot Gum')
figure(20202),plot([rhot,squeeze(KAPPAhatOOS(:,1,[4,6],1))]),legend('In sample GAS','OOS const','OOS GAS',2),title('Student''s t')
% very close. this is good - suggests that GAS parameters did not change much over the full sample relative to the IS

% now get the OOS copula log-likelihoods
outLLOOS = nan(P,6,2);
for uu=1:2;
    outLLOOS(:,1,uu) = NormalCopula_CLa(KAPPAhatIS(1,1,uu),UhatOOS(R+1:end,:,uu));
    outLLOOS(:,2,uu) = claytonCLa(KAPPAhatIS(1,2,uu),UhatOOS(R+1:end,:,uu));
    outLLOOS(:,3,uu) = gumbelCLa(KAPPAhatIS(1,3,uu),1-UhatOOS(R+1:end,:,uu));
    outLLOOS(:,4,uu) = tcopulaCL2a(KAPPAhatIS(1:2,4,uu),UhatOOS(R+1:end,:,uu));
    outLLOOS(:,5,uu) = Rotgumbel_GAS_CLa(KAPPAhatIS(1:3,5,uu),UhatOOS(R+1:end,:,uu),KAPPAhatIS(1,3,uu),[KK,HESSgumbel(:,1)]);
    outLLOOS(:,6,uu) = tcopula_GAS_CLa(KAPPAhatIS(1:4,6,uu),UhatOOS(R+1:end,:,uu),KAPPAhatIS(1,4,uu),RR,NN,HESSstudt);
end
outLLOOS = -outLLOOS;  % converting from neg log-like to regular log-like
sum(outLLOOS)

table8(7,:) = sum(outLLOOS(:,:,2));
table8(8,:) = ranks(-sum(outLLOOS(:,:,2))');
table8(8+7,:) = sum(outLLOOS(:,:,1));
table8(8+8,:) = ranks(-sum(outLLOOS(:,:,1))');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GW pair-wise tests for parametric and nonparametric margins

outLLGWtest = nan(6,6,2);
for uu=1:2
    for jj=1:5
        for ii=jj+1:6;
            temp = nwest(outLLOOS(:,ii,uu)-outLLOOS(:,jj,uu),ones(P,1),10);
            outLLGWtest(ii,jj,uu) = temp.tstat;
        end
    end
end
outLLGWtest
table8(1:6,1:6)  = outLLGWtest(:,:,2);
table8(9:14,1:6) = outLLGWtest(:,:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GW pair-wise test for mdoels with same copula but different margins (param vs nonparam)

% for EDF marginals I do not have a marginal density. For the density i need some sort of kernel (I guess). let's try using the built in one in
% Matlab for this 

mm=1;
[h1 f1 y1] = pltdens(stdresidsOOS(1:R,mm));
sum(f1<=0)  % 21
figure(10),plot(y1,f1,'b.-')

mm=2;
[h2 f2 y2] = pltdens(stdresidsOOS(1:R,mm));
sum(f2<=0)  % 29
figure(12),plot(y2,f2,'b.-')

edfpdf1 = interp1(y1,f1,stdresidsOOS(R+1:end,1));
min(edfpdf1)  %   0.00060259    so the neg value of the pdf turns out not to matter, so no problem
edfpdf2 = interp1(y2,f2,stdresidsOOS(R+1:end,2));
min(edfpdf2)  %   0.00033014    so the neg value of the pdf turns out not to matter, so no problem
% OK, so no problem with neg density. Go with this. The fact that the kernel density is negative (slightly) in some regions may cause problems in
% other applications so this needs to be monitored.

outLLcomp = nan(P,6+1,2);  % log-likelihood of joint models with nonparam margins vs param maregs (same copual)
outLLcomp(:,7,1) = log(edfpdf1)+log(edfpdf2);  % log like of the marginals
outLLcomp(:,7,2) = log(skewtdis_pdf(stdresidsOOS(R+1:end,1),outSKEWTIS(1,1),outSKEWTIS(1,2))) + log(skewtdis_pdf(stdresidsOOS(R+1:end,2),outSKEWTIS(2,1),outSKEWTIS(2,2)));
for uu=1:2;
    outLLcomp(:,1,uu) = -NormalCopula_CLa(KAPPAhatIS(1,1,uu),UhatOOS(R+1:end,:,uu));
    outLLcomp(:,2,uu) = -claytonCLa(KAPPAhatIS(1,2,uu),UhatOOS(R+1:end,:,uu));
    outLLcomp(:,3,uu) = -gumbelCLa(KAPPAhatIS(1,3,uu),1-UhatOOS(R+1:end,:,uu));
    outLLcomp(:,4,uu) = -tcopulaCL2a(KAPPAhatIS(1:2,4,uu),UhatOOS(R+1:end,:,uu));
    outLLcomp(:,5,uu) = -Rotgumbel_GAS_CLa(KAPPAhatIS(1:3,5,uu),UhatOOS(R+1:end,:,uu),KAPPAhatIS(1,3,uu),[KK,HESSgumbel(:,1)]);
    outLLcomp(:,6,uu) = -tcopula_GAS_CLa(KAPPAhatIS(1:4,6,uu),UhatOOS(R+1:end,:,uu),KAPPAhatIS(1,4,uu),RR,NN,HESSstudt);
end
squeeze(mean(outLLcomp))
squeeze(sum(outLLcomp))

outLLcompT = nan(1,6);
for cc=1:6;
    temp = nwest( outLLcomp(:,7,1)-outLLcomp(:,7,2)+outLLcomp(:,cc,1)-outLLcomp(:,cc,2),ones(P,1),10);
    outLLcompT(cc) = -temp.tstat;
end
table8(end,:) = outLLcompT;

info.fmt = strvcat('%10.2f');
info.cnames = strvcat('Normal','Clayton','Rot Gumbel','Stud t','RGum-GAS','Stud t-GAS');
info.rnames = strvcat('.','Normal','Clayton','Rot Gumbel','Stud t','RGum-GAS','Stud t-GAS','log L','Rank');
fprintf(['\n\nTable 8a: Out-of-sample model comparisons of copula models - parametric margins\n'])
mprint(table8(1:8,:),info)
fprintf(['\n\nTable 8b: Out-of-sample model comparisons of copula models - nonparametric margins\n'])
mprint(table8(9:16,:),info)
info.rnames = strvcat('.','t-stat');
fprintf(['\n\nTable 8c: Out-of-sample model comparisons of marginal models models - parametric vs nonparametric margins\n'])
mprint(table8(end,:),info)


% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_12.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_12.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR CORRELATION FROM COPULA-BASED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% let's first plot conditional tail depednence 
uu=1;
tauUg = zeros(T,1);
tauLg = 2-2.^(1./kappat);

tauUt = 2*tdis_cdf(-sqrt(1+1/theta8tvs(end))*sqrt((1-rhot)./(1+rhot)),1+1/theta8tvs(end));
tauLt = tauUt;

figure(1900),plot((1:T),tauLg,'r:','LineWidth',2);hold on;
plot((1:T),tauLt,'b-');
title('Tail dependence from time-varying copula models');
legend('RotGumbel lower tail','Stud t upper and lower tail',2);grid on;hold off;
set(gca,'XTick',jandates(1:2:end));
set(gca,'XTickLabel',datestr(datenum(datesYMD(jandates(1:2:end),1),datesYMD(jandates(1:2:end),2),datesYMD(jandates(1:2:end),3)),12));hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation approximation to linear correlation from a coupla-based models

% key here is for the estimate at each node to be VERY accurate. For a fixed computation time there is a trade off between the number of nodes and the
% number of simulations at each node. I suggest that 10-20 nodes with 100,000 simulations each works well (and better than 100 nodes with 10,000 simualtions each)

gridG1 = linspace(min(kappat), max(kappat), 10);
gridG2 = linspace(min(kappat), max(kappat), 20);

% WARNING: THIS PART IS PRETTY SLOW

sims = 1000;
sims2 = 100;  % doing many replications of length 1000 (could consider doing sims = 100,000 but may start to have memory problems)
outG1 = nan(length(gridG1),2,sims2);
tic;
for ss=1:sims2;
    for gg=1:length(gridG1)
        U3 = 1-Gumbel_rnd(gridG1(gg),sims);
        outG1(gg,1,ss) = corrcoef12(U3);  % rank correlation
        outG1(gg,2,ss) = corrcoef12([quantile(stdresids(:,1),U3(:,1),1),quantile(stdresids(:,2),U3(:,2),1)]);
    end
end
toc  % takes 12.2 mins for sims=1000 and sims2=100, length(gridG1)=10
outG1a = mean(outG1,3);

% WARNING: THIS PART IS PRETTY SLOW

outG2 = nan(length(gridG2),2);
tic;
for ss=1:sims2;
    for gg=1:length(gridG2)
        U3 = 1-Gumbel_rnd(gridG2(gg),sims);
        outG2(gg,1,ss) = corrcoef12(U3);  % rank correlation
        outG2(gg,2,ss) = corrcoef12([quantile(stdresids(:,1),U3(:,1),1),quantile(stdresids(:,2),U3(:,2),1)]);
    end
end
toc  % takes 24.4 mins for sims=1000 and sims2=100, length(gridG1)=20
outG2a = mean(outG2,3);

figure(1931),plot(gridG2,interp1(gridG1,outG1a(:,2),gridG2,'spline'),'-',gridG2,outG2a(:,2),'rp'),legend('Interpolated using 10 nodes','Values from 20 nodes',2),title('Linear correlation');
figure(1932),plot(gridG2,interp1(gridG1,outG1a(:,1),gridG2,'spline'),'-',gridG2,outG2a(:,1),'rp'),legend('Interpolated using 10 nodes','Values from 20 nodes',2),title('Rank correlation');


% now do same for Student's t
gridT1 = linspace(min(rhot), max(rhot), 10);
gridT2 = linspace(min(rhot), max(rhot), 20);

sims=100;
sim2=100;
outT1 = nan(length(gridT1),2,sims2);
outT2 = nan(length(gridT2),2,sims2);
tic;
for ss=1:sims2;
    for gg=1:length(gridT1)
        U4 = tdis_cdf(mvtrnd([[1,gridT1(gg)];[gridT1(gg),1]],1/theta8tvs(2),T),1/theta8tvs(2));
        outT1(gg,1,ss) = corrcoef12(U4);  % rank correlation
        outT1(gg,2,ss) = corrcoef12([quantile(stdresids(:,1),U4(:,1),1),quantile(stdresids(:,2),U4(:,2),1)]);
    end
end
toc  % takes 21 seconds for sims=100, sims2 = 100 and length(gridG1)=10
outT1a = mean(outT1,3);

tic;
for ss=1:sims2;
    for gg=1:length(gridT2)
        U4 = tdis_cdf(mvtrnd([[1,gridT2(gg)];[gridT2(gg),1]],1/theta8tvs(2),T),1/theta8tvs(2));
        outT2(gg,1,ss) = corrcoef12(U4);  % rank correlation
        outT2(gg,2,ss) = corrcoef12([quantile(stdresids(:,1),U4(:,1),1),quantile(stdresids(:,2),U4(:,2),1)]);
    end
end
toc  % takes 21 seconds for sims=100, sims2 = 100 and length(gridG1)=20
outT2a = mean(outT2,3);
  
figure(1951),plot(gridT2,interp1(gridT1,outT1a(:,2),gridT2,'spline'),'-',gridT2,outT2a(:,2),'rp'),legend('Interpolated using 10 nodes','Values from 20 nodes',2),title('Linear correlation');
figure(1952),plot(gridT2,interp1(gridT1,outT1a(:,1),gridT2,'spline'),'-',gridT2,outT2a(:,1),'rp'),legend('Interpolated using 10 nodes','Values from 20 nodes',2),title('Rank correlation');

figure(1953),subplot(2,1,1),plot(gridG2,interp1(gridG1,outG1a(:,2),gridG2,'spline'),'b-','LineWidth',2);hold on;
plot(gridG2,outG2a(:,2),'rp'),...
    title('Linear correlation from Gumbel copula model'),xlabel('Gumbel copula parameter'),ylabel('Linear correlation');grid on;hold off;
subplot(2,1,2),plot(gridT2,interp1(gridT1,outT1a(:,1),gridT2,'spline'),'b-','LineWidth',2);hold on;
plot(gridT2,outT2a(:,1),'rp'),legend('Interpolated using 10 nodes','Values from 20 nodes',2),...
    title('Linear correlation from Student''s t copula model'),xlabel('t copula correlation parameter'),ylabel('Linear correlation');grid on;hold off;

figure(1961),plot((1:T),interp1(gridG1,outG1a(:,2),kappat,'spline'),'r:','LineWidth',2);hold on;
plot((1:T),interp1(gridT1,outT1a(:,2),rhot,'spline'),'b-');
title('Linear correlation from time-varying copula models');
legend('RotGumbel','Stud t',2);grid on;hold off;
set(gca,'XTick',jandates(1:2:end));
set(gca,'XTickLabel',datestr(datenum(datesYMD(jandates(1:2:end),1),datesYMD(jandates(1:2:end),2),datesYMD(jandates(1:2:end),3)),12));hold off;
% looks good



% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_13.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_13.mat'';'];
evalin('base',temp_str);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VALUE-AT-RISK AND EXPECTED SHORTFALL FROM COPULA-BASED MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getting VaR and ES for a portfolio of these variables

% WARNING: THIS PART IS VERY SLOW

WW = [0;0.01;(0.05:0.05:0.95)';0.99;1];  % getting lots of weights, as this part is cheap once we have the siumulations
QQ = [0.001;0.005;0.01;0.05;0.1];  % quantiles to look at

muhat = rets1 - resids;  % conditional mean

sims = 5000;
outVAR = nan(T,2,length(WW),length(QQ),2,2);
tic;
for tt=1:T;
    U3 = 1-Gumbel_rnd(kappat(tt),sims);
    U4 = tdis_cdf(mvtrnd([[1,rhot(tt)];[rhot(tt),1]],1/theta8tvs(2),sims),1/theta8tvs(2));
    UUU = U3;
    UUU(:,:,2) = U4;
    EEE = nan(sims,2,2);
    YYY = nan(sims,2,2);
    for cc=1:2;
        for mm=1:2;
            EEE(:,mm,cc) = quantile(stdresids(:,mm),UUU(:,mm,cc),1);
            YYY(:,mm,cc) = muhat(tt,mm) + sqrt(hhat_opt(tt,mm))*EEE(:,mm,cc);  % simulated value for return
        end
        for ww=1:length(WW);
            w = WW(ww);
            pf = w*YYY(:,1,cc) + (1-w)*YYY(:,2,cc);
            pf2 = w*EEE(:,1,cc) + (1-w)*EEE(:,2,cc);
            outVAR(tt,cc,ww,:,1,1) = quantile(pf,QQ);                 % Value at Risk
            outVAR(tt,cc,ww,:,1,2) = quantile(pf2,QQ);                 % Value at Risk of the std resids (useful for seeing where the copula matters)
            
            for qq=1:length(QQ);
                temp124 = (pf<=quantile(pf,QQ(qq)));  % observations below this quantile
                if sum(temp124)>0
                    outVAR(tt,cc,ww,qq,2,1) = mean(pf(temp124));   % Expected Shortfall
                end
                temp124 = (pf2<=quantile(pf2,QQ(qq)));  % observations below this quantile
                if sum(temp124)>0
                    outVAR(tt,cc,ww,qq,2,2) = mean(pf2(temp124));   % Expected Shortfall
                end
            end
        end
    end
    if mod(tt,100)==0
        [tt,toc]
    end
end
toc  % takes about 3.88 hours for sims=5000

wstar = find(WW==0.5);
qq = find(QQ==0.01);
figure(2200+qq),subplot(2,1,1),plot((1:T)',squeeze(outVAR(:,1,wstar,qq,1,1)),'r:','LineWidth',2);hold on;
plot((1:T)',squeeze(outVAR(:,2,wstar,qq,1,1)),'b');
title(['Value-at-Risk from time-varying copula models, w=[',num2str(WW(wstar)),',',num2str(WW(wstar)),'], q=',num2str(QQ(qq))]);
legend('RotGumbel','Stud t',3);grid on;hold off;
set(gca,'XTick',jandates(1:2:end));
set(gca,'XTickLabel',datestr(datenum(datesYMD(jandates(1:2:end),1),datesYMD(jandates(1:2:end),2),datesYMD(jandates(1:2:end),3)),12));hold off;

subplot(2,1,2),plot((1:T)',squeeze(outVAR(:,1,wstar,qq,2,1)),'r:','LineWidth',2);hold on;
plot((1:T)',squeeze(outVAR(:,2,wstar,qq,2,1)),'b');
title(['Expected Shortfall from time-varying copula models, w=[',num2str(WW(wstar)),',',num2str(WW(wstar)),'], q=',num2str(QQ(qq))]);
legend('RotGumbel','Stud t',3);grid on;hold off;
set(gca,'XTick',jandates(1:2:end));
set(gca,'XTickLabel',datestr(datenum(datesYMD(jandates(1:2:end),1),datesYMD(jandates(1:2:end),2),datesYMD(jandates(1:2:end),3)),12));hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% let's look at VAR and ES as a function of rank correlation for these two copulas (cleaner mapping, and since there's only one parametric varying in
% each of these models this will make a nice clean figure 

% WARNING: THIS PART IS PRETTY SLOW

gridG3 = [1;1.05;1.1;linspace(1.2, max(kappat)*1.5, 17)'];  % adding some more poitns near zero as there is a lot of curvature there

sims = 1000;
sims2 = 100;  % doing many replications of length 1000 (could consider doing sims = 100,000 but may start to have memory problems)
outG3 = nan(length(gridG3),2,sims2);
tic;
for ss=1:sims2;
    for gg=1:length(gridG3)
        U3 = 1-Gumbel_rnd(gridG3(gg),sims);
        outG3(gg,1,ss) = corrcoef12(U3);  % rank correlation
        outG3(gg,2,ss) = corrcoef12([quantile(stdresids(:,1),U3(:,1),1),quantile(stdresids(:,2),U3(:,2),1)]);
    end
end
toc  % takes 12.2 mins for sims=1000 and sims2=100, length(gridG1)=10
outG3a = mean(outG3,3);
figure(934234),plot(gridG3,outG3a,'bo-')

gridT3 = linspace(0, 0.99, 20);
sims=100;
sim2=100;
outT3 = nan(length(gridT3),2,sims2);
tic;
for ss=1:sims2;
    for gg=1:length(gridT3)
        U4 = tdis_cdf(mvtrnd([[1,gridT3(gg)];[gridT3(gg),1]],1/theta8tvs(2),T),1/theta8tvs(2));
        outT3(gg,1,ss) = corrcoef12(U4);  % rank correlation
        outT3(gg,2,ss) = corrcoef12([quantile(stdresids(:,1),U4(:,1),1),quantile(stdresids(:,2),U4(:,2),1)]);
    end
end
toc  % takes 21 seconds for sims=100, sims2 = 100 and length(gridG1)=10
outT3a = mean(outT3,3);
figure(9234),plot(gridT3,outT3a,'bo-')

RR = [(0:0.05:0.95)';0.99];  % covering the range of rank corrlations 
WW = [0;0.01;(0.05:0.05:0.95)';0.99;1];  % getting lots of weights, as this part is cheap once we have the siumulations
QQ = [0.001;0.005;0.01;0.05;0.1];  % quantiles to look at

sims = 5000;

outVAR2 = nan(length(RR),3,length(WW),length(QQ),2,sims2);  % rank correl ; [gumbel, t copula, Normal] ; pf weight ; quantile ; [VaR, ES]

sims = 5000;
sims2 = 500;
tic;
for ss=1:sims2;
    for rr=1:length(RR);
        gstar = interp1(outG3a(:,1),gridG3',RR(rr),'spline');  % value of Gumbel model that best matches this rank correlation, be inverting using mapping established above
        rstar = interp1(outT3a(:,1),gridT3',RR(rr),'spline');
        U3 = 1-Gumbel_rnd(gstar,sims);
        U4 = tdis_cdf(mvtrnd([[1,rstar];[rstar,1]],1/theta8tvs(2),sims),1/theta8tvs(2));
        U5 = normcdf(mvnrnd(zeros(1,2),[[1,rstar];[rstar,1]],sims));
        UUU = U3;
        UUU(:,:,2) = U4;
        UUU(:,:,3) = U5;
        EEE = nan(sims,2,3);
        for cc=1:3; 
            for mm=1:2;
                EEE(:,mm,cc) = quantile(stdresids(:,mm),UUU(:,mm,cc),1);
            end
            for ww=1:length(WW);
                w = WW(ww);
                pf2 = w*EEE(:,1,cc) + (1-w)*EEE(:,2,cc);
                outVAR2(rr,cc,ww,:,1,ss) = quantile(pf2,QQ);                 % Value at Risk of the std resids (useful for seeing where the copula matters)
                
                for qq=1:length(QQ);
                    temp124 = (pf2<=quantile(pf2,QQ(qq)));  % observations below this quantile
                    if sum(temp124)>0
                        outVAR2(rr,cc,ww,qq,2,ss) = mean(pf2(temp124));   % Expected Shortfall
                    end
                end
            end
        end
        [ss,rr,toc]  % takes 34 seconds per loop for sims=50,000. Takes ?? seconds for sims=100,000
    end
end  
outVAR2a = mean(outVAR2,6);
toc  % takes 1.05 hours for sims=5000, sims2=100

wstar = find(WW==0.5);
qq = find(QQ==0.01);
figure(2660+qq),subplot(2,2,1),plot(RR,squeeze(outVAR2a(:,3,wstar,qq,1,1)),'kp-');hold on;
plot(RR,squeeze(outVAR2a(:,2,wstar,qq,1,1)),'b.-');
plot(RR,squeeze(outVAR2a(:,1,wstar,qq,1)),'ro:','LineWidth',2);
%axis([0,1,-4.3,-1.7]);
axis([0,1,-2.65,-1.75]);
title(['Portfolio Value-at-Risk, q=',num2str(QQ(qq))]);
xlabel('Rank correlation');
grid on;hold off;
subplot(2,2,3),plot(RR,squeeze(outVAR2a(:,3,wstar,qq,2)),'kp-','LineWidth',1);hold on;
plot(RR,squeeze(outVAR2a(:,2,wstar,qq,2,1)),'b.-');
plot(RR,squeeze(outVAR2a(:,1,wstar,qq,2)),'ro:','LineWidth',2);hold on;
%axis([0,1,-5.1,-2.1]);
axis([0,1,-3.3,-2.1]);
title(['Portfolio Expected Shortfall, q=',num2str(QQ(qq))]);
xlabel('Rank correlation');
grid on;hold off;
qq = find(QQ==0.001);
subplot(2,2,2),plot(RR,squeeze(outVAR2a(:,3,wstar,qq,1,1)),'kp-');hold on;
plot(RR,squeeze(outVAR2a(:,2,wstar,qq,1,1)),'b.-');
plot(RR,squeeze(outVAR2a(:,1,wstar,qq,1)),'ro:','LineWidth',2);
%axis([0,1,-4.3,-1.7]);
axis([0,1,-4.3,-2.6]);
title(['Portfolio Value-at-Risk, q=',num2str(QQ(qq))]);
xlabel('Rank correlation');
legend('Normal','Stud t','RotGumbel',1);grid on;hold off;
subplot(2,2,4),plot(RR,squeeze(outVAR2a(:,3,wstar,qq,2)),'kp-','LineWidth',1);hold on;
plot(RR,squeeze(outVAR2a(:,2,wstar,qq,2,1)),'b.-');
plot(RR,squeeze(outVAR2a(:,1,wstar,qq,2)),'ro:','LineWidth',2);hold on;
%axis([0,1,-5.1,-2.1]);
axis([0,1,-5.1,-3]);
title(['Portfolio Expected Shortfall, q=',num2str(QQ(qq))]);
xlabel('Rank correlation');
grid on;hold off;



% saving the output down to this part of the code.
temp_str = ['save ''',save_path,save_name,'_stage_14.mat'';'];
evalin('base',temp_str);
% to load the output matrix from this point run the following
temp_str = ['load ''',save_path,save_name,'_stage_14.mat'';'];
evalin('base',temp_str);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%