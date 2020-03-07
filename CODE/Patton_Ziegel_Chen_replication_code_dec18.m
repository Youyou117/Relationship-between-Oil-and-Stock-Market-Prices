% Code to replicate the results for the S&P 500 index returns in: 
%
% Patton, A.J., J.F. Ziegel, and R. Chen, 2018, Dynamic Semiparametric
% Models for Expected Shortfall (and Value-at-Risk), Journal of Econometrics, forthcoming. 
%
% 17 December 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loading in the data

input_file_path = 'D:\Work\Johanna_Rui\';  % change this to the directory where you save the data file

data1 = load([input_file_path,'sp500_rets_1990_2016.txt'],'ascii');
data1 = data1((data1(:,2)~=0) ,:);    % dropping zero returns (where market was closed)
% columns are date (YYYYMMDD) and daily return (in percent)
T = length(data1)  % 6800 obs

split_date = 20000000;  % use 20000000 to make the last day of the in-sample period 19991231.
data2 = data1(data1(:,1)<split_date,2);  % only using data from 1990-1999 (first ten years of sample)
R = length(data2);
P = T-R;
data3 = data1(data1(:,1)>split_date,2);  % out-of-sample period is all obs after 19991231.

figure(1),plot(data1(:,2))

options = optimset('Display','iter','MaxFunEvals',2000);  % changing default options so that I can see iterations, and increasing number of function evals (fminsearch can require many of these)

VEhatALL = nan(P,3+3+4,2);  % models are RW-125, RW-250, RW-500, GARCH-N, GARCH-skew t, GARCH-EDF, FZ2F, FZ1F, GARCH-FZ, Hybrid

alpha = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% summary stats (Table 6)
AA = [0.01;0.025;0.05;0.10];
summ_stats = nan(15,1);
summ_stats(1:3) = [min(data1(:,1)); max(data1(:,1)); length(data1) ];
summ_stats(4:7) = [252*mean(data1(:,2)),sqrt(252)*std(data1(:,2)),skewness(data1(:,2)),kurtosis(data1(:,2))];
for aa=1:length(AA)
    summ_stats([7+aa,7+aa+length(AA)],1) = sample_VE(data1(:,2),AA(aa));
end
summ_stats

clear info;
info.cnames = strvcat('SP500');
info.rnames = strvcat('.','Mean','StdDev','Skew','Kurt','VaR-0.01','VaR-0.025','VaR-0.05','VaR-0.10','ES-0.01','ES-0.025','ES-0.05','ES-0.10');
info.fmt    = '%10.3f';
mprint(summ_stats(4:end,:),info)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mean, variance, skew t estimates (generates the entries in Table 7)

% find the optimal ARMA model
maxP=6;
maxQ=6;

tic;
[theta,sig2,vcv,order,resids] = ARMAX_opt(data2,maxP,maxQ,'BIC');
orderALL = order;
toc  % takes 11 seconds
orderALL  % optimal order is (0,0), i.e., just a constant

phi0 = mean(data2)

% estimate the GARCH(1,1) parameters
lower = [2.1,-0.99]';
upper = [50,0.99]';
theta0 = [5,-0.1]';
GARCHparamALL = nan(3+2,1);  % will save the estimated parameters.

tic;
% estimating a GARCH model on the residuals
resids2 = data2-mean(data2);
[parameters, ~, hhat2 ] = garchpq(resids2, 1, 1);  % need Kevin Sheppard's GARCH toolbox for this step
GARCHparamALL(1:3,1) = parameters';

% estimating a skew t distribution on the std resids
skewtparams = fmincon('skewtdis_LL',theta0,[],[],[],[],lower,upper,[],options,resids2./sqrt(hhat2));
GARCHparamALL(4:5,1) = skewtparams;
toc  % 1.3 seconds
GARCHparamALL

% getting OOS forecasts for GARCH models
muhatOOS = nan(P,1);
muhatOOS = mean(data2)*ones(P,1);  % easy: for this asset the conditional mean was estimated as a constant

hhatOOS = nan(P,1);
hhatOOS(1) = GARCHparamALL(1:3,1)'*[1;(data2(end)-mean(data2))^2;hhat2(end)];  % first obs of OOS period relies on last obs of in-sample period
for tt=2:P
    hhatOOS(tt) = GARCHparamALL(1:3,1)'*[1;(data3(tt-1)-muhatOOS(tt-1))^2;hhatOOS(tt-1)];
end

VEhatALL(:,4,:) = muhatOOS + sqrt(hhatOOS)*[norminv(alpha),-normpdf(norminv(alpha))/alpha];  % VaR and ES forecasts assuming std Normal innovations
VEhatALL(:,5,:) = muhatOOS + sqrt(hhatOOS)*skewt_VE(GARCHparamALL(4:5),alpha);               % VaR and ES forecasts assuming skew t innovations
VEhatALL(:,6,:) = muhatOOS + sqrt(hhatOOS)*sample_VE(resids2./sqrt(hhat2),alpha);            % VaR and ES forecasts using the empirical dist of std resids

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getting OOS forecasts using rolling windows

MM = [125;250;500];
tic;
for mm=1:length(MM)
    M = MM(mm);
    for tt=R+1:T
        VEhatALL(tt-R,mm,:) = sample_VE(data1(tt-M:tt-1,2),alpha);
    end
end
toc  % 0.2 seconds


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimating GAS models, in-sample (generates the entries in Table 8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GARCH-FZ
garchFZparam = nan(4,2);  % parameter estimates and std errors
tic;
% first starting value: the QLIKE GARCH parameters
VEbar = sample_VE( resids2./sqrt(hhat2) , alpha)';
theta0 = [norminv(GARCHparamALL(3));  log(GARCHparamALL(2)); log(-VEbar(2)) ; norminv(VEbar(1)/VEbar(2)) ];  % transforming the parameters to match the GARCH-FZ loss function.
[thetahat1,LL1] = fminsearch('garch_FZ_LL',theta0,options,resids2,alpha,[],GARCHparamALL(1));  % using the omega from QLIKE estimation to fix omega here
[theta0, thetahat1]

% second starting value: use smoothed value of FZ loss function
warning off;
[thetahat2,LL2] = fminunc('garch_FZ_LL',   theta0,   options,resids2,alpha,5,GARCHparamALL(1));
[thetahat3,LL3] = fminunc('garch_FZ_LL',   thetahat2,options,resids2,alpha,20,GARCHparamALL(1));
[thetahat4,LL4] = fminsearch('garch_FZ_LL',thetahat3,options,resids2,alpha,[],GARCHparamALL(1));
[theta0,thetahat2,thetahat3,thetahat4]
[LL1,LL4]
if LL1<LL4  % then using QLIKE results as starting values is preferred
    thetahat5 = thetahat1;
else
    thetahat5 = thetahat4;
end
garchFZparam = thetahat5;

[Eloss,VEhat] = feval('garch_FZ_LL',thetahat5,resids2,alpha,[],GARCHparamALL(1));
figure(2),plot(VEhat(:,3))

% now getting standard errors and "original" (untransformed) parameters
[ElossGFZ,VEhatGFZ,~,VCVgarchFZ,garchFZparamT] = garch_FZ_LL(garchFZparam,resids2,alpha,[],GARCHparamALL(1),[],length(resids2)^(-1/3));
garchFZparamT
ElossGFZ
toc % takes 1.94 seconds

% getting OOS forecasts
[~,fhat] = garch_FZ_LL(garchFZparam,data1(:,2)-mean(data2),alpha,[],GARCHparamALL(1),cov(resids2)); % getting fitted VaR and ES, using in-sample parameter estimate and starting value
VEhatALL(:,9,:) = muhatOOS + fhat(R+1:end,1:2);  % OOS VaR and ES forecasts from this model, adding back the mean

figure(1),plot(squeeze(VEhatALL(:,:,2)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one-factor GAS model

FZ1Fparam = nan(4,1);

tic;
VEbar = sample_VE( data2 , alpha)';
theta0 = [norminv(0.95);  log(0.005); log(-VEbar(2)) ; norminv(VEbar(1)/VEbar(2)) ];  % transforming the parameters to match the GARCH-FZ loss function
warning off;
[thetahat2,LL2] = fminunc('GAS_onefactor_LL3',theta0,       options,data2,alpha,5);
[thetahat3,LL3] = fminunc('GAS_onefactor_LL3',   thetahat2, options,data2,alpha,20);
[thetahat4,LL4] = fminsearch('GAS_onefactor_LL3',thetahat3, options,data2,alpha,[]);
FZ1Fparam = thetahat4;

[ElossFZ1F,VEhatFZ1F,loss,VCVFZ1F,FZ1FparamT] = GAS_onefactor_LL3(FZ1Fparam,data2,alpha,[],length(data2)^(-1/3));VCVFZ1F,FZ1FparamT(:,1),FZ1FparamT
ElossFZ1F
toc

% getting OOS forecasts
[~,fhat] = GAS_onefactor_LL3(FZ1Fparam,data1(:,2),alpha); % getting fitted VaR and ES, using in-sample parameter estimate and starting value
VEhatALL(:,8,:) = fhat(R+1:end,1:2);  % OOS VaR and ES forecasts from this model, adding back the mean

figure(1),plot(squeeze(VEhatALL(:,[1,6,8,9],2)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hybrid model

tic;
theta0 = [norminv(0.95);  log(0.005);   log(0.005); log(-VEbar(2)) ; norminv(VEbar(1)/VEbar(2)) ];  % transforming the parameters to match the GARCH-FZ loss function
warning off;
[thetahat2,LL2] = fminunc(   'GAS_hybrid_LL3',theta0,   options,data2,alpha,5);
[thetahat3,LL3] = fminunc(   'GAS_hybrid_LL3',thetahat2,options,data2,alpha,20);
[thetahat4,LL4] = fminsearch('GAS_hybrid_LL3',thetahat3,options,data2,alpha,[]);
[thetahat5,LL5] = fminsearch('GAS_hybrid_LL3',[FZ1Fparam(1:2);-10;FZ1Fparam(3:4)],options,data2,alpha,[]);          % using the FZ1F model as the starting value.
if LL4<LL5
    hybridFZparam = thetahat4;
else
    hybridFZparam = thetahat5;
end
    
% now getting standard errors and "original" (untransformed) parameters
[ElossHYBRID,VEhatHYBRID,~,VCVhybrid,hybridFZparamT] = GAS_hybrid_LL3(hybridFZparam,data2,alpha,[],length(data2)^(-1/3));
hybridFZparamT
ElossHYBRID
toc % takes 1.34 seconds


[~,fhat] = GAS_hybrid_LL3(hybridFZparam,data1(:,2),alpha,[],[],mean(log(abs(data2)))); % getting fitted VaR and ES, using in-sample parameter estimate and starting value
VEhatALL(:,10,:) = fhat(R+1:end,1:2);  % OOS VaR and ES forecasts from this model, adding back the mean

figure(1),plot(squeeze(VEhatALL(:,[1,6,8,9,10],2)));legend('RW-125','G-EDF','FZ1F','GFZ','Hybrid');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% two-factor GAS model
tic;
options = optimset('Display','off','MaxFunEvals',2000);  % changing default options so that I can see iterations, and increasing number of function evals (fminsearch can require many of these)

% this is a big model, so I'm trying a few different starting values
theta00 = [...
    [ (1-0.97)*VEbar(1) ; (1-0.97)*VEbar(2) ; norminv(0.97) ; norminv(0.97) ; 0      ; 0      ; 0      ;  0     ] , ...
    [ (1-0.98)*VEbar(1) ; (1-0.98)*VEbar(2) ; norminv(0.98) ; norminv(0.98) ; 0      ; 0      ; 0      ;  0     ] , ...
    [ (1-0.99)*VEbar(1) ; (1-0.99)*VEbar(2) ; norminv(0.99) ; norminv(0.99) ; 0      ; 0      ; 0      ;  0     ] ];

warning off;
outSTART = nan(size(theta00,1)+1,size(theta00,2));
for ss = 1:size(theta00,2)
    theta0 = theta00(:,ss);
    
    temp = GAS_twofactor_LL3(theta0,data2,alpha,[],5);
    if temp<1e6 && ~isnan(temp) && ~isinf(temp)  % then it's an OK starting value
        [thetahat2,LL2] = fminunc('GAS_twofactor_LL3',theta0,options,data2,alpha,[],5);
    else
        thetahat2 = theta0;
    end
    
    temp = GAS_twofactor_LL3(thetahat2,data2,alpha,[],20);
    if temp<1e6 && ~isnan(temp) && ~isinf(temp)  % then it's an OK starting value
        thetahat2a = thetahat2;
    else
        thetahat2a = theta0;  % reset to the original starting value
    end
    [thetahat3,LL3] = fminunc('GAS_twofactor_LL3',thetahat2a,   options,data2,alpha,[],20);
    
    temp = GAS_twofactor_LL3(thetahat3,data2,alpha,[],[]);
    if temp<1e6 && ~isnan(temp) && ~isinf(temp)  % then it's an OK starting value
        thetahat3a = thetahat3;
    else
        thetahat3a = theta0;  % reset to the original starting value
    end
    [thetahat4,LL4] = fminsearch('GAS_twofactor_LL3',thetahat3a,options,data2,alpha,[],[]);
    
    outSTART(:,ss) = [thetahat4;LL4];
    [ss,toc]  % takes about 25 seconds per starting value
end
outSTART
factor2FZparam = outSTART(1:end-1,find(outSTART(end,:)==min(outSTART(end,:)),1));  % picking the best parameter
min(outSTART(end,:))
[Eloss2FZ,VEhat2FZ,~,~,VCV2FZ,factor2FZparamT] = GAS_twofactor_LL3(factor2FZparam,data2,alpha,[],[],length(data2)^(-1/3));
factor2FZparamT
toc  % 45 seconds


[~,fhat] = GAS_twofactor_LL3(factor2FZparam,data1(:,2),alpha); % getting fitted VaR and ES, using in-sample parameter estimate and starting value
VEhatALL(:,7,:) = fhat(R+1:end,1:2);  % OOS VaR and ES forecasts from this model, adding back the mean

figure(1),plot(squeeze(VEhatALL(:,[1,6,7,8,10],2)))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
% OOS comparison of the VaR and ES forecasts

squeeze(sum(isnan(VEhatALL)))'
[VEhatALLa,out2] = VEforecast_filter(VEhatALL);
nansum(out2)  % no insane forecasts. Great!

% Diebold-Mariano tests
outDM = nan(10+1,10);
for ii=1:10
    [Eloss1, loss1] = VaR_ES_loss_0(squeeze(VEhatALLa(:,ii,:)),data3,alpha);
    outDM(end,ii) = Eloss1;
    for jj=ii+1:10
        [Eloss2, loss2] = VaR_ES_loss_0(squeeze(VEhatALLa(:,jj,:)),data3,alpha);
        temp = nwest(loss1-loss2,ones(P,1),20);
        outDM(ii,jj) = temp.tstat;
        outDM(jj,ii) = -temp.tstat;
    end
end


% Mincer-Zarnowitz (type) tests
lamVsALL = nan(P,10);
lamEsALL = nan(P,10);
outMZ = nan(10,2);
for mm=1:10
    lamVsALL(:,mm) = (data3<VEhatALLa(:,mm,1)) - alpha;
    lamEsALL(:,mm) = 1/alpha*(data3<VEhatALLa(:,mm,1)).*data3./VEhatALLa(:,mm,2) - 1;
    
    temp = nwest(lamVsALL(2:end,mm),[ones(P-1,1),lamVsALL(1:end-1,mm),VEhatALLa(2:end,mm,1)],20);
    outMZ(mm,1) = 1-chi2cdf(temp.beta'*inv(temp.vcv)*temp.beta,length(temp.beta));
    
    temp = nwest(lamEsALL(2:end,mm),[ones(P-1,1),lamEsALL(1:end-1,mm),VEhatALLa(2:end,mm,2)],20);
    outMZ(mm,2) = 1-chi2cdf(temp.beta'*inv(temp.vcv)*temp.beta,length(temp.beta));
end


clear info;
info.cnames = strvcat('Avg loss','MZ-VaR','MZ-ES');
info.rnames = strvcat('.','RW-125','RW-250','RW-500','GCH-n','GCH-skt','GCH-edf','FZ2F','FZ1F','GCH-FZ','Hybrid');
info.fmt    = '%10.3f';
info.width = 250;
sprintf('\n\n Table 9\n\n')
mprint([outDM(end,:)',outMZ(:,1),outMZ(:,2)],info)

info.cnames = strvcat('RW-125','RW-250','RW-500','GCH-n','GCH-skt','GCH-edf','FZ2F','FZ1F','GCH-FZ','Hybrid');
info.rnames = strvcat('.','RW-125','RW-250','RW-500','GCH-n','GCH-skt','GCH-edf','FZ2F','FZ1F','GCH-FZ','Hybrid');
sprintf('\n\n Table 10\n\n')
mprint(outDM(1:end-1,:),info)

  

%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% some figures

% full-sample VaR and ES: RW-125, GARCH-EDF, FZ1F
VEhatISOOS = nan(T,3,2);
M = 125;
for tt=M+1:T
    VEhatISOOS(tt,1,:) = sample_VE(data1(tt-M:tt-1,2),alpha);
end
VEhatISOOS(:,2,:) = mean(data2) + sqrt([hhat2(1:end);hhatOOS])*sample_VE(resids2./sqrt(hhat2),alpha);
[~,fhat] = GAS_onefactor_LL3(FZ1Fparam,data1(:,2),alpha); % getting fitted VaR and ES, using in-sample parameter estimate and starting value
VEhatISOOS(:,3,:) = fhat(:,1:2);

figure;plot(VEhatISOOS(:,:,2))

dates = data1(:,1);
dates2 = datesYMD(dates);
years = (1990:2016)';
jandates = 1;
for yy=2:length(years)
    temp = min(find(dates2(:,1)==years(yy)));
    jandates = [jandates;temp];
end
jandates = [jandates;T];


ymin1 = -8;
ymin2 = -12;
figure(101)
subplot(2,1,1),plot((1:T)',VEhatISOOS(:,3,1),'b-','LineWidth',1);hold on;
plot((1:T)',VEhatISOOS(:,2,1),'r--','LineWidth',1);
plot((1:T)',VEhatISOOS(:,1,1),'Color',0.65*ones(1,3),'LineWidth',2);
legend('One-factor GAS','GARCH-EDF','RW-125','Location','SouthWest');grid on;
xlim([0,T+1]);ylim([ymin1,0]);
title('5% Value-at-Risk forecasts for S&P 500 daily returns');
set(gca,'XTick',jandates([(1:3:end-2),end]));  
set(gca,'XTickLabel',datestr(datenum(dates2(jandates([(1:3:end-2),end]),1),dates2(jandates([(1:3:end-2),end]),2),dates2(jandates([(1:3:end-2),end]),3)),12));

subplot(2,1,2),plot((1:T)',VEhatISOOS(:,3,2),'b-','LineWidth',1);hold on;
plot((1:T)',VEhatISOOS(:,2,2),'r--','LineWidth',1);
plot((1:T)',VEhatISOOS(:,1,2),'Color',0.65*ones(1,3),'LineWidth',2);grid on;
xlim([0,T+1]);ylim([ymin2,0]);
title('5% Expected Shortfall forecasts for S&P 500 daily returns');
set(gca,'XTick',jandates([(1:3:end-2),end]));  
set(gca,'XTickLabel',datestr(datenum(dates2(jandates([(1:3:end-2),end]),1),dates2(jandates([(1:3:end-2),end]),2),dates2(jandates([(1:3:end-2),end]),3)),12));



jandates2 = [];
for jj=1:12
    jandates2 = [jandates2;min(find( (dates2(:,1)==2015).*(dates2(:,2)==jj) ))];
    jandates2 = [jandates2;min(find( (dates2(:,1)==2016).*(dates2(:,2)==jj) ))];
end
jandates2 = sort(jandates2);
jandates2 = [jandates2;length(dates2)];


temp124 = find(dates>20150000);
ymin1 = -4;
ymin2 = -5;
figure(102)
subplot(2,1,1),plot(temp124,VEhatISOOS(temp124,3,1),'b-','LineWidth',1);hold on;
plot(temp124,VEhatISOOS(temp124,2,1),'r--','LineWidth',1);
plot(temp124,VEhatISOOS(temp124,1,1),'Color',0.65*ones(1,3),'LineWidth',2);
legend('One-factor GAS','GARCH-EDF','RW-125','Location','SouthWest');grid on;
xlim([min(temp124)-1,max(temp124)+1]);ylim([ymin1,0]);
title('5% Value-at-Risk forecasts for S&P 500 daily returns');
set(gca,'XTick',jandates2([(1:3:end-2),end]));  
set(gca,'XTickLabel',datestr(datenum(dates2(jandates2([(1:3:end-2),end]),1),dates2(jandates2([(1:3:end-2),end]),2),dates2(jandates2([(1:3:end-2),end]),3)),12));

subplot(2,1,2),plot(temp124,VEhatISOOS(temp124,3,2),'b-','LineWidth',1);hold on;
plot(temp124,VEhatISOOS(temp124,2,2),'r--','LineWidth',1);
plot(temp124,VEhatISOOS(temp124,1,2),'Color',0.65*ones(1,3),'LineWidth',2);grid on;
xlim([min(temp124)-1,max(temp124)+1]);ylim([ymin2,0]);
title('5% Expected Shortfall forecasts for S&P 500 daily returns');
set(gca,'XTick',jandates2([(1:3:end-2),end]));  
set(gca,'XTickLabel',datestr(datenum(dates2(jandates2([(1:3:end-2),end]),1),dates2(jandates2([(1:3:end-2),end]),2),dates2(jandates2([(1:3:end-2),end]),3)),12));
