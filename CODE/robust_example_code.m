% This program replications the empirical results and figures using the "robust" family
% of loss functions for volatility forecast comparison proposed in:
%
% Patton, A.J., 2006, "Volatility Forecast Comparison using Imperfect Volatility Proxies", 
% manuscript, London School of Economics.
%
%  Andrew Patton
%
%  11 August, 2007.

% running this code took 28.5 seconds on my 2.66GHz machine with 3GB of RAM
%
% (it only takes a few seconds for the first two parts: the little
% simulation in section 3 below takes the most time)

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 1.  PLOTTING THE LOSS FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r2=2;
inc = 0.01;
hhat = (inc:inc:4)';
hhat1 = (2-inc:-inc:inc)';
hhat2 = (2+inc:inc:4-inc)';
hhat3 = (inc:inc:2-inc)';
figure(1),plot(...
    hhat,robust_loss_1(1,r2,hhat),':',...
    hhat,robust_loss_1(0.5,r2,hhat),'.-',...
    hhat,robust_loss_1(0,r2,hhat),...
    hhat,robust_loss_1(-0.5,r2,hhat),'-.',...
    hhat,robust_loss_1(-1,r2,hhat),...
    hhat,robust_loss_1(-2,r2,hhat),'--',...
    hhat,robust_loss_1(-5,r2,hhat),'-x',...
    [r2,r2],[0,100],'k--'),...
    axis([0,4,0,2.5]),...
    legend('b=1','b=0.5','b=0 (MSE)','b=-0.5','b=-1','b=-2 (QLIKE)','b=-5'),...
    xlabel('hhat (r2=2)'),ylabel('loss'),...
    title('Various robust loss functions');

figure(2),plot(...
    hhat3,robust_loss_1(1,r2,hhat2)./robust_loss_1(1,r2,hhat1),':',...
    hhat3,robust_loss_1(0.5,r2,hhat2)./robust_loss_1(0.5,r2,hhat1),'.-',...
    [0,4],[1,1],...
    hhat3,robust_loss_1(-0.5,r2,hhat2)./robust_loss_1(-0.5,r2,hhat1),'-.',...
    hhat3,robust_loss_1(-1,r2,hhat2)./robust_loss_1(-1,r2,hhat1),...
    hhat3,robust_loss_1(-2,r2,hhat2)./robust_loss_1(-2,r2,hhat1),'--',...
    hhat3,robust_loss_1(-5,r2,hhat2)./robust_loss_1(-5,r2,hhat1),'-x',...
    [r2,r2],[0,100],'k--'),...
    axis([0,2,0,2.5]),...
    legend('b=1','b=0.5','b=0 (MSE)','b=-0.5','b=-1','b=-2 (QLIKE)','b=-5'),...
    xlabel('forecast error (r2=2)'),ylabel('loss'),...
    title('Ratio of loss from negative forecast errors to positive forecast errors');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 2.  EMPIRICAL RESULTS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load 'robust_ibm_data_apr06.txt' -ascii;

dates = robust_ibm_data_apr06(:,1);     % dates in YYYY.MMDD format. from 1993.0104 to 2003.1231
rets = robust_ibm_data_apr06(:,2);      % daily returns
RV = robust_ibm_data_apr06(:,3:end);    % realised vol(m=1,m=6, m=26, m=78 corresponding to daily, 65-min, 15-min, 5-min samping)

T = size(rets,1);
R = length(find(dates<19940000));  % "in-sample" period = 253  obs
n = T-R; % out-of-sample period 2519 obs
n=2500   % increasing the in-sample by 19 days so that the out-of-sample is exactly 2500 obs
R = T-n

lam = 0.94;
ww = 60;  % window length for the rolling window estimator
HHAT = nines(T,2);  % hhat1 = riskmetrics , hhat2 rolling window 
HHAT(1:R,1) = mean(RV(1:R,1))*ones(R,1);  % mean squared daily return over in-sample period ~~ sample variance that year
HHAT(1:R,2) = mean(RV(R-(ww-1):R,1))*ones(R,1);  % average 5-min RV over the last WW(ww) days

for tt=R+1:T;
    HHAT(tt,1) = 0.94*HHAT(tt-1,1) + 0.06*RV(tt-1,1);
    HHAT(tt,2) = mean(RV(tt-ww:tt-1,1));
end


% Mincer-Zarnowitz regressions
outMZ = nines(2+2+2+1,size(HHAT,2),size(RV,2));  % b0hat, b1hat, b0se, b1se, chi2stat, pval, R2
for ii=1:size(HHAT,2);
    for jj=1:size(RV,2);
        temp = nwest(RV(R+1:end,jj),[ones(n,1),HHAT(R+1:end,ii)]);
        chi2stat = (temp.beta-[0;1])'*inv(temp.vcv)*(temp.beta-[0;1]);
        pval = 1-chi2cdf(chi2stat,2);
        outMZ(:,ii,jj) = [temp.beta;temp.se;chi2stat;pval;temp.rsqr];
    end
end
% in the format of Table 3 in the paper:
format bank;
squeeze(outMZ([1;3;2;4;5;6],2,:))   % results for rolling window
squeeze(outMZ([1;3;2;4;5;6],1,:))   % results for RiskMetrics


% now doing some Diebold-Mariano-West tests
BB = [1;0;-1;-2;-5];    %loss function shape parameters
outDMW = nines(length(BB),size(RV,2));
for jj=1:size(RV,2);
    [out1 out2] = robust_loss_1(BB,RV(:,jj),HHAT);
    outDMW(:,jj) = squeeze(out2(2,1,:));  % want to have rolling window as hhat1, so take the (2,1) element of out2
end
% Table 4 in the paper:
outDMW
    

%%%%% plotting the forecasts

dates2a = dates2(dates(R+1:end)*10000);		% transforming the date format
% getting the row numbers corresponding to first day of each year in sample
jandates = nines(11,3);
index1 = nines(11,1);
for ii = 1:10;
   index1(ii) = min(find((dates2a(:,2)==1).*(dates2a(:,1)==1993+ii)));
   jandates(ii,:) = dates2a(index1(ii),:);
end
index1(end) = n;
jandates(end,:) = dates2a(end,:);

figure(3),plot((1:n)',RV(R+1:end,1),'g-',(1:n)',RV(R+1:end,end),'k-',...
    (1:n)',HHAT(R+1:end,2),'r-',(1:n),HHAT(R+1:end,1),'b-'),...
    title('Conditional variance forecasts'),ylabel('Conditional variance'),...
    legend('daily squared returns','RV-5min','60-day rolling window','RiskMetrics',2),...
    set(gca,'XTick',index1),...
    set(gca,'XTickLabel',datestr(datenum(dates2a(index1,1),dates2a(index1,2),dates2a(index1,3)),12));
% too messy!

figure(4),plot((1:n)',HHAT(R+1:end,2),'r-','LineWidth',2');hold on;
plot((1:n),HHAT(R+1:end,1),'b-'),...
    title('Conditional variance forecasts'),ylabel('Conditional variance'),...
    legend('60-day rolling window','RiskMetrics',2),...
    set(gca,'XTick',index1),...
    set(gca,'XTickLabel',datestr(datenum(dates2a(index1,1),dates2a(index1,2),dates2a(index1,3)),12));hold off;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% 3.  THE DMW STATISTIC EXAMPLE    %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% here I present some a simulation study to support the simple DMW 
% exmample in the paper

a = 0.05;  
b = 0.90;
w = 1-a-b;
% can alter the parameters above (note that w does not matter so long as it
% is positive) but need the following to be positive, for E[sig[t]^4] to exist 
1 - (b^2) - 2*a*b - 3*(a^2)  

NN = [50;100;150;250;500;1000;2500;5000;10000];  % some sample sizes
reps = 100;  

out1 = nines(reps,length(NN),2);
for ii=1:reps;
    for jj=1:length(NN);
        n = NN(jj);
        [simulatedata, H] = garchsimulate(n,[w;a;b],1,1);
        hhat1 = H;
        hhat2 = 2/pi*H;
        volproxy = simulatedata.^2;
        dt = ((sqrt(volproxy)-sqrt(hhat1)).^2) - ((sqrt(volproxy)-sqrt(hhat2)).^2);
        temp = ols(dt,ones(n,1));  % naive standard errors
        out1(ii,jj,1) = temp.tstat;
        temp = nwest(dt,ones(n,1),ceil(n^(1/3)));  % robust standard errors
        out1(ii,jj,2) = temp.tstat;
    end
    if mod(ii,10)==0
        [ii,toc/60]
    end
end
toc
NN1 = (min(NN):1:max(NN))';
figure(5),plot(NN,squeeze(mean(out1(:,:,1))),'o-',...
    NN,squeeze(mean(out1(:,:,2))),'v-');hold on;
plot(NN1,1/sqrt((5+3*sqrt(2/pi))/(1-sqrt(2/pi))*(1-((a+b)^2))/(1-((a+b)^2)-2*(a^2)) - 1)*sqrt(NN1),'r:','LineWidth',2);
plot(NN,1.96*ones(length(NN),1),'k--'),...
    legend('naive std errors','robust std errors','theoretical',2),...
    xlabel('sample size'),ylabel('DMW_0 statistic'),...
    title('Diebold-Mariano-West statistic for MSD-SD loss example');hold off;