% Replicating the results in:
%
% Bollerslev, T., A. J. Patton, and R. Quaedvlieg, 2016, 
% Exploiting the Errors: A Simple Approach for Improved Volatility Forecasting, 
% Journal of Econometrics, 192, 1-18.
%
%  Andrew Patton and Rogier Quaedvlieg
%
%  20 December 2016

% NOTE: The computations for the published paper were all done in Ox. The
% code below produces numbers that differ slightly from those in the
% published article. We provide Matlab code for greater accessibility.
%
% In addition, some of the results in the published article were
% inaccurate. The corrected tables are in an accompanying PDF file. The
% interpretation of the results remains valid, even after the correcting
% for the innacuracies.  


% This code takes about 3.5 mins to run completely. Almost all of this is
% in the rolling and expanding window estimations, each of which takes
% about 10-15 seconds to run.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first: load in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

temp = xlsread('E:\BPQ2016_Data\Market\SP500_RV_5min.xlsx','A2:L4097');  % change this line and the one below to where you have saved the spreadsheets of data
% Tables 9 and 10 use data on the SPY, not the SP500 index like the rest of the tables.
data1min = xlsread('E:\BPQ2016_Data\Market\SPY_RV_1min.xlsx','1','b2:g4203');  % 1 RV, 2 SSRV, 3 TSRV, 4 Kernel, 5 PARV, 6 RQ


dates = temp(:,1);
data = temp(:,2:end);  % RV	RQ	RJ	BPV	RV_min	RV_plus	TPQ	MedRQ	TrRQ	RQ15min	 RQboot

RV  = data(:,1);
RQ  = data(:,2);
RJ  = data(:,3);
BPV = data(:,4);
RVn = data(:,5);
RVp = data(:,6);

TPQ     = data(:,7);
MedRQ   = data(:,8);
TrRQ    = data(:,9);
RQ15min = data(:,10);
RQboot  = data(:,11);

clear temp;
T = length(RV);

insanityfilter = 1;  % =1 if want filter, equal 0 if not


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 2: summary statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table2 = nan(1,6); 
table2(1:4) = [min(RV),mean(RV),median(RV),max(RV)];
table2(5) = corrcoef12(RV(2:end),RV(1:end-1));

% ARQ
sqrtRQdemean = sqrt(RQ) - mean(sqrt(RQ));  % we use this in our "HARQ" model
tempARQ = hwhite(RV(23:end),[ones(T-22,1),RV(22:end-1),RV(22:end-1).*sqrtRQdemean(22:end-1)]);
table2(6) = tempARQ.beta(2);

clear info;
info.rnames = char('Table 2','SP500');
info.cnames = char('Min','Mean','Median','Max','AR','ARQ');
info.fmt    = '%10.3f';
sprintf('\n\n');
mprint(table2,info)

% NOTE: we used Hansen-Lunde's Ox code for the "AR-HL" estimator. See their
% web page for details and code. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 3: in-sample results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RVd = RV;
RVw = mean([RV,mlag(RV,4,mean(RV))],2);
RVm = mean([RV,mlag(RV,21,mean(RV))],2);  % note these are not lagged, so need to lag one period if using in a forecasting regression

sqrtRQdemeanD = sqrt(RQ) - mean(sqrt(RQ));  % we use this in our "HARQ" model
temp = mean([RQ,mlag(RQ,4,mean(RQ))],2);
sqrtRQdemeanW = sqrt(temp) - mean(sqrt(temp));
temp = mean([RQ,mlag(RQ,21,mean(RQ))],2);
sqrtRQdemeanM = sqrt(temp) - mean(sqrt(temp));

% AR(1)
tempAR = hwhite(RV(23:end),[ones(T-22,1),RV(22:end-1)]);
hhatAR = [ones(T-22,1),RV(22:end-1)]*tempAR.beta;  % dropping first 21 obs to make comparable with HAR

% HAR
tempHAR = hwhite(RV(23:end),[ones(T-22,1),RVd(22:end-1),RVw(22:end-1),RVm(22:end-1)]);
hhatHAR = [ones(T-22,1),RVd(22:end-1),RVw(22:end-1),RVm(22:end-1)]*tempHAR.beta;

% HARQ
tempHARQ = hwhite(RV(23:end),[ones(T-22,1),RVd(22:end-1),RVw(22:end-1),RVm(22:end-1),RVd(22:end-1).*sqrtRQdemeanD(22:end-1)]);
hhatHARQ = [ones(T-22,1),RVd(22:end-1),RVw(22:end-1),RVm(22:end-1),RVd(22:end-1).*sqrtRQdemeanD(22:end-1)]*tempHARQ.beta;

% ARQ
tempARQ = hwhite(RV(23:end),[ones(T-22,1),RVd(22:end-1),RVd(22:end-1).*sqrtRQdemeanD(22:end-1)]);
hhatARQ = [ones(T-22,1),RVd(22:end-1),RVd(22:end-1).*sqrtRQdemeanD(22:end-1)]*tempARQ.beta;

% HARQF
tempHARQF = hwhite(RV(23:end),[ones(T-22,1),RVd(22:end-1),RVw(22:end-1),RVm(22:end-1),RVd(22:end-1).*sqrtRQdemeanD(22:end-1),RVw(22:end-1).*sqrtRQdemeanW(22:end-1),RVm(22:end-1).*sqrtRQdemeanM(22:end-1)]);
hhatHARQF = [ones(T-22,1),RVd(22:end-1),RVw(22:end-1),RVm(22:end-1),RVd(22:end-1).*sqrtRQdemeanD(22:end-1),RVw(22:end-1).*sqrtRQdemeanW(22:end-1),RVm(22:end-1).*sqrtRQdemeanM(22:end-1)]*tempHARQF.beta;


table3 = nan(7*2+3,5);
table3([1,3],1) = tempAR.beta;
table3([2,4],1) = tempAR.se;
table3((1:2:7),2) = tempHAR.beta;
table3((2:2:8),2) = tempHAR.se;
table3([1,3,9],3) = tempARQ.beta;
table3([2,4,10],3) = tempARQ.se;
table3((1:2:9),4) = tempHARQ.beta;
table3((2:2:10),4) = tempHARQ.se;
table3((1:2:13),5) = tempHARQF.beta;
table3((2:2:14),5) = tempHARQF.se;

hhatALLis = [hhatAR,hhatHAR,hhatARQ,hhatHARQ,hhatHARQF];
for ii=1:size(hhatALLis,2);
    table3(15,ii) = 1-cov(RV(23:end)-hhatALLis(:,ii))/cov(RV(23:end));            % r2
    table3(16,ii) = mean( (RV(23:end)-hhatALLis(:,ii)).^2);                              % MSE
    table3(17,ii) = mean( RV(23:end)./hhatALLis(:,ii) - log(RV(23:end)./hhatALLis(:,ii)) -1 );     % QLIKE
end

clear info;
info.cnames = char('AR','HAR','ARQ','HARQ','HARQF');
info.rnames = char('Table 3','b0','s.e.','b1','s.e.','b2','s.e.','b3','s.e.','b1Q','s.e.','b2Q','s.e.','b3Q','s.e.','R2','MSE','QLIKE');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table3,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 4: out-of-sample results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, the rolling window results
window = 1000;
yhatALLrw = nan(T,1+8);  % will store the realized value and the forecasts: AR, HAR, HAR-J, CHAR, SHAR, ARQ, HARQ, HARQF
countrw = zeros(8,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(tt-1000:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(tt-1000:tt-1);
    
    minRV = min(RV1);  % will use these values in the "insanity filter" below
    maxRV = max(RV1);
    meanRV = mean(RV1);
 
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    RV1 = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtRQw = sqrt(mean([RQ1,mlag(RQ1, 4,mean(RQ1))],2));
    sqrtRQm = sqrt(mean([RQ1,mlag(RQ1,21,mean(RQ1))],2));
    sqrtRQd = sqrtRQd(22:end);
    sqrtRQw = sqrtRQw(22:end);
    sqrtRQm = sqrtRQm(22:end);
    
    RJ1  = RJ(tt-1000+21:tt-1);
    RVp1 = RVp(tt-1000+21:tt-1);
    RVn1 = RVn(tt-1000+21:tt-1);

    BPV1  = BPV(tt-1000:tt-1);
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1,4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);

    % getting the realised vol for this day
    yhatALLrw(tt,1) = RV(tt);
    
    % 1 AR
    tempAR = ols(RV1(2:end),[ones(window-22,1),RV1(1:end-1)]);
    yhatALLrw(tt,1+1) = [1,RV1(end)]*tempAR.beta;
    
    % 2 HAR
    tempHAR = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHAR.beta;

    % 3 HAR-J
    tempHARJ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RJ1(1:end-1)]);
    yhatALLrw(tt,1+3) = [1,RV1d(end),RV1w(end),RV1m(end),RJ1(end)]*tempHARJ.beta;
    
    % 4 CHAR
    tempCHAR = ols(RV1(2:end),[ones(window-22,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1)]);
    yhatALLrw(tt,1+4) = [1,BPV1d(end),BPV1w(end),BPV1m(end)]*tempCHAR.beta;
    
    % 5 SHAR
    tempSHAR = ols(RV1(2:end),[ones(window-22,1),RVp1(1:end-1),RVn1(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLrw(tt,1+5) = [1,RVp1(end),RVn1(end),RV1w(end),RV1m(end)]*tempSHAR.beta;
    
    % 6 ARQ
    tempARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+6) = [1,RV1d(end),RV1d(end).*sqrtRQd(end)]*tempARQ.beta;
    
    % 7 HARQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+7) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 8 HARQF
    tempHARQF = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1),RV1w(1:end-1).*sqrtRQw(1:end-1),RV1m(1:end-1).*sqrtRQm(1:end-1)]);
    yhatALLrw(tt,1+8) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end),RV1w(end).*sqrtRQw(end),RV1m(end).*sqrtRQm(end)]*tempHARQF.beta;
        
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:8
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+8);  % will store the realized value and the forecasts
countiw = zeros(8,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(1:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(1:tt-1);
    
    minRV = min(RV(tt-1000:tt-1));  % will use these values in the "insanity filter" below
    maxRV = max(RV(tt-1000:tt-1));
    meanRV = mean(RV(tt-1000:tt-1));
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    RV1  = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtRQw = sqrt(mean([RQ1,mlag(RQ1, 4,mean(RQ1))],2));
    sqrtRQm = sqrt(mean([RQ1,mlag(RQ1,21,mean(RQ1))],2));
    sqrtRQd = sqrtRQd(22:end);
    sqrtRQw = sqrtRQw(22:end);
    sqrtRQm = sqrtRQm(22:end);
    
    RJ1  = RJ(1+21:tt-1);
    RVp1 = RVp(1+21:tt-1);
    RVn1 = RVn(1+21:tt-1);

    BPV1  = BPV(1:tt-1);
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1,4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);
    
    % getting the realised vol for this day
    yhatALLiw(tt,1) = RV(tt);
    
    % 1 AR
    tempAR = ols(RV1(2:end),[ones(tt-23,1),RV1(1:end-1)]);
    yhatALLiw(tt,1+1) = [1,RV1(end)]*tempAR.beta;
    
    % 2 HAR
    tempHAR = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHAR.beta;

    % 3 HAR-J
    tempHARJ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RJ1(1:end-1)]);
    yhatALLiw(tt,1+3) = [1,RV1d(end),RV1w(end),RV1m(end),RJ1(end)]*tempHARJ.beta;
    
    % 4 CHAR
    tempCHAR = ols(RV1(2:end),[ones(tt-23,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1)]);
    yhatALLiw(tt,1+4) = [1,BPV1d(end),BPV1w(end),BPV1m(end)]*tempCHAR.beta;
    
    % 5 SHAR
    tempSHAR = ols(RV1(2:end),[ones(tt-23,1),RVp1(1:end-1),RVn1(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLiw(tt,1+5) = [1,RVp1(end),RVn1(end),RV1w(end),RV1m(end)]*tempSHAR.beta;
    
    % 6 ARQ
    tempARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+6) = [1,RV1d(end),RV1d(end).*sqrtRQd(end)]*tempARQ.beta;
    
    % 7 HARQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+7) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 8 HARQF
    tempHARQF = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1),RV1w(1:end-1).*sqrtRQw(1:end-1),RV1m(1:end-1).*sqrtRQm(1:end-1)]);
    yhatALLiw(tt,1+8) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end),RV1w(end).*sqrtRQw(end),RV1m(end).*sqrtRQm(end)]*tempHARQF.beta;
        
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:8
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countiw(ii) = countiw(ii) + temp;
        end
    end
end


yhatALL = nan(T,1+8,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1000 obs where we have no OOS forecast

table4a = nan(4,8);
table4 = nan(4,8);
for mm=[2,1,(3:8)]
    table4a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table4a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table4a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table4a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window

    table4(:,mm) = table4a(:,mm)./table4a(:,2);
end
clear info;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F');
info.rnames = char('Table 4','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table4,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 5: stratified OOS losses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RQ1 = RQ(window:end-1); % pulling out RQ for the OOS period.
table5a = nan(4+4,8);
table5 = nan(4+4,8);
for mm=[2,1,(3:8)]
    temp124 = find(RQ1<=quantile(RQ1,0.95)); % days with low RQ
    temp125 = find(RQ1>quantile(RQ1,0.95));  % days with high RQ
    
    table5a(1,mm) = mean( (yhatALL(temp124,1,1)-yhatALL(temp124,mm+1,1)).^2 );  % MSE - rolling window
    table5a(2,mm) = mean( (yhatALL(temp124,1,2)-yhatALL(temp124,mm+1,2)).^2 );  % MSE - expanding window
    table5a(3,mm) = mean( yhatALL(temp124,1,1)./yhatALL(temp124,mm+1,1) - log(yhatALL(temp124,1,1)./yhatALL(temp124,mm+1,1)) -1 );  % QLIKE - rolling window
    table5a(4,mm) = mean( yhatALL(temp124,1,2)./yhatALL(temp124,mm+1,2) - log(yhatALL(temp124,1,2)./yhatALL(temp124,mm+1,2)) -1 );  % QLIKE - expanding window

    table5a(5,mm) = mean( (yhatALL(temp125,1,1)-yhatALL(temp125,mm+1,1)).^2 );  % MSE - rolling window
    table5a(6,mm) = mean( (yhatALL(temp125,1,2)-yhatALL(temp125,mm+1,2)).^2 );  % MSE - expanding window
    table5a(7,mm) = mean( yhatALL(temp125,1,1)./yhatALL(temp125,mm+1,1) - log(yhatALL(temp125,1,1)./yhatALL(temp125,mm+1,1)) -1 );  % QLIKE - rolling window
    table5a(8,mm) = mean( yhatALL(temp125,1,2)./yhatALL(temp125,mm+1,2) - log(yhatALL(temp125,1,2)./yhatALL(temp125,mm+1,2)) -1 );  % QLIKE - expanding window

    table5(:,mm) = table5a(:,mm)./table5a(:,2);
end
clear info;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F');
info.rnames = char('.','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\nTable 5: Bottom 95% RQ');
mprint(table5(1:4,:),info)
sprintf('Table 5: Top 5% RQ');
mprint(table5(5:8,:),info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 6: In-sample weekly and monthly estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table6 = nan(2*7,4+4);
T = size(RV,1);
RVd = RV;
RVw = mean([RV,mlag(RV,4,mean(RV))],2);
RVm = mean([RV,mlag(RV,21,mean(RV))],2);  % note these are not lagged, so need to lag one period if using in a forecasting regression

sqrtRQdemeanD = sqrt(RQ) - mean(sqrt(RQ));  % we use this in our "HARQ" model
temp = mean([RQ,mlag(RQ,4,mean(RQ))],2); sqrtRQdemeanW = sqrt(temp) - mean(sqrt(temp));
temp = mean([RQ,mlag(RQ,21,mean(RQ))],2);sqrtRQdemeanM = sqrt(temp) - mean(sqrt(temp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h=5

% HAR
tempHAR = nwest(RVw(23+4:end),[ones(T-22-4,1),RVd(22:end-1-4),RVw(22:end-1-4),RVm(22:end-1-4)],10);

% HARQ
tempHARQ = nwest(RVw(23+4:end),[ones(T-22-4,1),RVd(22:end-1-4),RVw(22:end-1-4),RVm(22:end-1-4),RVd(22:end-1-4).*sqrtRQdemeanD(22:end-1-4)],10);

% HARQF
tempHARQF = nwest(RVw(23+4:end),[ones(T-22-4,1),RVd(22:end-1-4),RVw(22:end-1-4),RVm(22:end-1-4),RVd(22:end-1-4).*sqrtRQdemeanD(22:end-1-4),RVw(22:end-1-4).*sqrtRQdemeanW(22:end-1-4),RVm(22:end-1-4).*sqrtRQdemeanM(22:end-1-4)],10);

% HARQF-m
tempHARQh = nwest(RVw(23+4:end),[ones(T-22-4,1),RVd(22:end-1-4),RVw(22:end-1-4),RVm(22:end-1-4),RVw(22:end-1-4).*sqrtRQdemeanW(22:end-1-4)],10);

table6(1:2:8,1) = tempHAR.beta;
table6(2:2:8,1) = tempHAR.se;
table6(1:2:10,2) = tempHARQ.beta;
table6(2:2:10,2) = tempHARQ.se;
table6(1:2:14,3) = tempHARQF.beta;
table6(2:2:14,3) = tempHARQF.se;
table6([(1:2:8)';11],4) = tempHARQh.beta;
table6([(2:2:8)';12],4) = tempHARQh.se;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h=22

% HAR
tempHAR = nwest(RVm(23+21:end),[ones(T-22-21,1),RVd(22:end-1-21),RVw(22:end-1-21),RVm(22:end-1-21)],44);

% HARQ
tempHARQ = nwest(RVm(23+21:end),[ones(T-22-21,1),RVd(22:end-1-21),RVw(22:end-1-21),RVm(22:end-1-21),RVd(22:end-1-21).*sqrtRQdemeanD(22:end-1-21)],44);

% HARQF
tempHARQF = nwest(RVm(23+21:end),[ones(T-22-21,1),RVd(22:end-1-21),RVw(22:end-1-21),RVm(22:end-1-21),RVd(22:end-1-21).*sqrtRQdemeanD(22:end-1-21),RVw(22:end-1-21).*sqrtRQdemeanW(22:end-1-21),RVm(22:end-1-21).*sqrtRQdemeanM(22:end-1-21)],44);

% HARQF-m
tempHARQh = nwest(RVm(23+21:end),[ones(T-22-21,1),RVd(22:end-1-21),RVw(22:end-1-21),RVm(22:end-1-21),RVm(22:end-1-21).*sqrtRQdemeanM(22:end-1-21)],44);

table6(1:2:8,5) = tempHAR.beta;
table6(2:2:8,5) = tempHAR.se;
table6(1:2:10,6) = tempHARQ.beta;
table6(2:2:10,6) = tempHARQ.se;
table6(1:2:14,7) = tempHARQF.beta;
table6(2:2:14,7) = tempHARQF.se;
table6([(1:2:8)';13],8) = tempHARQh.beta;
table6([(2:2:8)';14],8) = tempHARQh.se;

clear info;
info.cnames = char('HAR','HARQ','HARQ-F','HARQ-h','HAR','HARQ','HARQ-F','HARQ-h');
info.rnames = char('Table 6','b0','s.e.','b2','s.e.','b3','s.e.','b4','s.e.','b1Q','s.e.','b2Q','s.e.','b3Q','s.e.');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table6,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 7: OOS weekly estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, the rolling window results
window = 1000;
yhatALLrw = nan(T,1+9);  % will store the realized value and the forecasts
countrw = zeros(9,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(tt-1000:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(tt-1000:tt-1);

    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1, 4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    minRV = min(RV1w);  % will use these values in the "insanity filter" below
    maxRV = max(RV1w);
    meanRV = mean(RV1w);

    RV1  = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtRQw = sqrt(mean([RQ1,mlag(RQ1, 4,mean(RQ1))],2));
    sqrtRQm = sqrt(mean([RQ1,mlag(RQ1,21,mean(RQ1))],2));
    sqrtRQd = sqrtRQd(22:end);
    sqrtRQw = sqrtRQw(22:end);
    sqrtRQm = sqrtRQm(22:end);
    
    RJ1  =  RJ(tt-1000+21:tt-1);
    RVp1 = RVp(tt-1000+21:tt-1);
    RVn1 = RVn(tt-1000+21:tt-1);

    BPV1  = BPV(tt-1000:tt-1);
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1, 4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);
    
    % getting the realised vol for this day
    yhatALLrw(tt,1) = RVw(tt);
    
    % 1 AR
    tempAR = ols(RV1w(4+2:end),[ones(window-22-4,1),RV1d(1:end-1-4)]);
    yhatALLrw(tt,1+1) = [1,RV1d(end-4)]*tempAR.beta;
    
    % 2 HAR
    tempHAR = ols(RV1w(4+2:end),[ones(window-22-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4)]*tempHAR.beta;

    % 3 HAR-J
    tempHARJ = ols(RV1w(4+2:end),[ones(window-22-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RJ1(1:end-1-4)]);
    yhatALLrw(tt,1+3) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RJ1(end-4)]*tempHARJ.beta;
    
    % 4 CHAR
    tempCHAR = ols(RV1w(4+2:end),[ones(window-22-4,1),BPV1d(1:end-1-4),BPV1w(1:end-1-4),BPV1m(1:end-1-4)]);
    yhatALLrw(tt,1+4) = [1,BPV1d(end-4),BPV1w(end-4),BPV1m(end-4)]*tempCHAR.beta;
    
    % 5 SHAR
    tempSHAR = ols(RV1w(4+2:end),[ones(window-22-4,1),RVp1(1:end-1-4),RVn1(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4)]);
    yhatALLrw(tt,1+5) = [1,RVp1(end-4),RVn1(end-4),RV1w(end-4),RV1m(end-4)]*tempSHAR.beta;
    
    % 6 ARQ
    tempARQ = ols(RV1w(4+2:end),[ones(window-22-4,1),RV1d(1:end-1-4),RV1d(1:end-1-4).*sqrtRQd(1:end-1-4)]);
    yhatALLrw(tt,1+6) = [1,RV1d(end-4),RV1d(end-4).*sqrtRQd(end-4)]*tempARQ.beta;
    
    % 7 HARQ
    tempHARQ = ols(RV1w(4+2:end),[ones(window-22-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RV1d(1:end-1-4).*sqrtRQd(1:end-1-4)]);
    yhatALLrw(tt,1+7) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RV1d(end-4).*sqrtRQd(end-4)]*tempHARQ.beta;
    
    % 8 HARQF
    tempHARQF = ols(RV1w(4+2:end),[ones(window-22-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RV1d(1:end-1-4).*sqrtRQd(1:end-1-4),RV1w(1:end-1-4).*sqrtRQw(1:end-1-4),RV1m(1:end-1-4).*sqrtRQm(1:end-1-4)]);
    yhatALLrw(tt,1+8) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RV1d(end-4).*sqrtRQd(end-4),RV1w(end-4).*sqrtRQw(end-4),RV1m(end-4).*sqrtRQm(end-4)]*tempHARQF.beta;

    % 9 HARQh
    tempHARQh = ols(RV1w(4+2:end),[ones(window-22-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RV1w(1:end-1-4).*sqrtRQw(1:end-1-4)]);
    yhatALLrw(tt,1+9) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RV1w(end-4).*sqrtRQw(end-4)]*tempHARQh.beta;
        
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:9
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+9);  % will store the realized value and the forecasts
countiw = zeros(9,1);
for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(1:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(1:tt-1);
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    minRV = min(RV1w);  % will use these values in the "insanity filter" below
    maxRV = max(RV1w);
    meanRV = mean(RV1w);
    
    RV1  = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtRQw = sqrt(mean([RQ1,mlag(RQ1, 4,mean(RQ1))],2));
    sqrtRQm = sqrt(mean([RQ1,mlag(RQ1,21,mean(RQ1))],2));
    sqrtRQd = sqrtRQd(22:end);
    sqrtRQw = sqrtRQw(22:end);
    sqrtRQm = sqrtRQm(22:end);
    
    RJ1  = RJ(1+21:tt-1);
    RVp1 = RVp(1+21:tt-1);
    RVn1 = RVn(1+21:tt-1);

    BPV1  = BPV(1:tt-1);
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1,4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);
    
    % getting the realised vol for this day
    yhatALLiw(tt,1) = RVw(tt);
    
    % 1 AR
    tempAR = ols(RV1w(4+2:end),[ones(tt-23-4,1),RV1(1:end-1-4)]);
    yhatALLiw(tt,1+1) = [1,RV1(end-4)]*tempAR.beta;
    
    % 2 HAR
    tempHAR = ols(RV1w(4+2:end),[ones(tt-23-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4)]*tempHAR.beta;

    % 3 HAR-J
    tempHARJ = ols(RV1w(4+2:end),[ones(tt-23-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RJ1(1:end-1-4)]);
    yhatALLiw(tt,1+3) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RJ1(end-4)]*tempHARJ.beta;
    
    % 4 CHAR
    tempCHAR = ols(RV1w(4+2:end),[ones(tt-23-4,1),BPV1d(1:end-1-4),BPV1w(1:end-1-4),BPV1m(1:end-1-4)]);
    yhatALLiw(tt,1+4) = [1,BPV1d(end-4),BPV1w(end-4),BPV1m(end-4)]*tempCHAR.beta;
    
    % 5 SHAR
    tempSHAR = ols(RV1w(4+2:end),[ones(tt-23-4,1),RVp1(1:end-1-4),RVn1(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4)]);
    yhatALLiw(tt,1+5) = [1,RVp1(end-4),RVn1(end-4),RV1w(end-4),RV1m(end-4)]*tempSHAR.beta;
    
    % 6 ARQ
    tempARQ = ols(RV1w(4+2:end),[ones(tt-23-4,1),RV1d(1:end-1-4),RV1d(1:end-1-4).*sqrtRQd(1:end-1-4)]);
    yhatALLiw(tt,1+6) = [1,RV1d(end-4),RV1d(end-4).*sqrtRQd(end-4)]*tempARQ.beta;
    
    % 7 HARQ
    tempHARQ = ols(RV1w(4+2:end),[ones(tt-23-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RV1d(1:end-1-4).*sqrtRQd(1:end-1-4)]);
    yhatALLiw(tt,1+7) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RV1d(end-4).*sqrtRQd(end-4)]*tempHARQ.beta;
    
    % 8 HARQF
    tempHARQF = ols(RV1w(4+2:end),[ones(tt-23-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RV1d(1:end-1-4).*sqrtRQd(1:end-1-4),RV1w(1:end-1-4).*sqrtRQw(1:end-1-4),RV1m(1:end-1-4).*sqrtRQm(1:end-1-4)]);
    yhatALLiw(tt,1+8) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RV1d(end-4).*sqrtRQd(end-4),RV1w(end-4).*sqrtRQw(end-4),RV1m(end-4).*sqrtRQm(end-4)]*tempHARQF.beta;
        
    % 9 HARQF
    tempHARQh = ols(RV1w(4+2:end),[ones(tt-23-4,1),RV1d(1:end-1-4),RV1w(1:end-1-4),RV1m(1:end-1-4),RV1w(1:end-1-4).*sqrtRQw(1:end-1-4)]);
    yhatALLiw(tt,1+9) = [1,RV1d(end-4),RV1w(end-4),RV1m(end-4),RV1w(end-4).*sqrtRQw(end-4)]*tempHARQh.beta;
 
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:9
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countiw(ii) = countiw(ii) + temp;
        end
    end
end 

yhatALL = nan(T,1+9,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1001 obs where we have no OOS forecast

table7a = nan(4,9);
table7 = nan(4,9);
for mm=[2,1,(3:9)]
    table7a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table7a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table7a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table7a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window

    table7(:,mm) = table7a(:,mm)./table7a(:,2);
end
clear info;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F','HARQ-h');
info.rnames = char('Table 7','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
info.width = 150;
sprintf('\n\n');
mprint(table7,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 8: OOS monthly estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, the rolling window results
window = 1000;
yhatALLrw = nan(T,1+9);  % will store the realized value and the forecasts
countrw = zeros(9,1);


for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(tt-1000:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(tt-1000:tt-1);

    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1, 4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    minRV = min(RV1m);  % will use these values in the "insanity filter" below
    maxRV = max(RV1m);
    meanRV = mean(RV1m);

    RV1  = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtRQw = sqrt(mean([RQ1,mlag(RQ1, 4,mean(RQ1))],2));
    sqrtRQm = sqrt(mean([RQ1,mlag(RQ1,21,mean(RQ1))],2));
    sqrtRQd = sqrtRQd(22:end);
    sqrtRQw = sqrtRQw(22:end);
    sqrtRQm = sqrtRQm(22:end);
    
    RJ1  =  RJ(tt-1000+21:tt-1);
    RVp1 = RVp(tt-1000+21:tt-1);
    RVn1 = RVn(tt-1000+21:tt-1);

    BPV1  = BPV(tt-1000:tt-1);
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1, 4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);
    
    % getting the realised vol for this day
    yhatALLrw(tt,1) = RVm(tt);
    
    % 1 AR
    tempAR = ols(RV1m(21+2:end),[ones(window-22-21,1),RV1d(1:end-1-21)]);
    yhatALLrw(tt,1+1) = [1,RV1d(end-21)]*tempAR.beta;
    
    % 2 HAR
    tempHAR = ols(RV1m(21+2:end),[ones(window-22-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21)]*tempHAR.beta;

    % 3 HAR-J
    tempHARJ = ols(RV1m(21+2:end),[ones(window-22-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RJ1(1:end-1-21)]);
    yhatALLrw(tt,1+3) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RJ1(end-21)]*tempHARJ.beta;
    
    % 4 CHAR
    tempCHAR = ols(RV1m(21+2:end),[ones(window-22-21,1),BPV1d(1:end-1-21),BPV1w(1:end-1-21),BPV1m(1:end-1-21)]);
    yhatALLrw(tt,1+4) = [1,BPV1d(end-21),BPV1w(end-21),BPV1m(end-21)]*tempCHAR.beta;
    
    % 5 SHAR
    tempSHAR = ols(RV1m(21+2:end),[ones(window-22-21,1),RVp1(1:end-1-21),RVn1(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21)]);
    yhatALLrw(tt,1+5) = [1,RVp1(end-21),RVn1(end-21),RV1w(end-21),RV1m(end-21)]*tempSHAR.beta;
    
    % 6 ARQ
    tempARQ = ols(RV1m(21+2:end),[ones(window-22-21,1),RV1d(1:end-1-21),RV1d(1:end-1-21).*sqrtRQd(1:end-1-21)]);
    yhatALLrw(tt,1+6) = [1,RV1d(end-21),RV1d(end-21).*sqrtRQd(end-21)]*tempARQ.beta;
    
    % 7 HARQ
    tempHARQ = ols(RV1m(21+2:end),[ones(window-22-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RV1d(1:end-1-21).*sqrtRQd(1:end-1-21)]);
    yhatALLrw(tt,1+7) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RV1d(end-21).*sqrtRQd(end-21)]*tempHARQ.beta;
    
    % 8 HARQF
    tempHARQF = ols(RV1m(21+2:end),[ones(window-22-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RV1d(1:end-1-21).*sqrtRQd(1:end-1-21),RV1w(1:end-1-21).*sqrtRQw(1:end-1-21),RV1m(1:end-1-21).*sqrtRQm(1:end-1-21)]);
    yhatALLrw(tt,1+8) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RV1d(end-21).*sqrtRQd(end-21),RV1w(end-21).*sqrtRQw(end-21),RV1m(end-21).*sqrtRQm(end-21)]*tempHARQF.beta;

    % 9 HARQh
    tempHARQh = ols(RV1m(21+2:end),[ones(window-22-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RV1m(1:end-1-21).*sqrtRQm(1:end-1-21)]);
    yhatALLrw(tt,1+9) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RV1m(end-21).*sqrtRQm(end-21)]*tempHARQh.beta;
        
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:9
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+9);  % will store the realized value and the forecasts
countiw = zeros(9,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(1:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(1:tt-1);
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    minRV = min(RV1m);  % will use these values in the "insanity filter" below
    maxRV = max(RV1m);
    meanRV = mean(RV1m);
    
    RV1 = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtRQw = sqrt(mean([RQ1,mlag(RQ1, 4,mean(RQ1))],2));
    sqrtRQm = sqrt(mean([RQ1,mlag(RQ1,21,mean(RQ1))],2));
    sqrtRQd = sqrtRQd(22:end);
    sqrtRQw = sqrtRQw(22:end);
    sqrtRQm = sqrtRQm(22:end);
    
    RJ1  = RJ(1+21:tt-1);
    RVp1 = RVp(1+21:tt-1);
    RVn1 = RVn(1+21:tt-1);
    
    BPV1  = BPV(1:tt-1);
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1,4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);

    % getting the realised vol for this day
    yhatALLiw(tt,1) = RVm(tt);
    
    % 1 AR
    tempAR = ols(RV1m(21+2:end),[ones(tt-23-21,1),RV1(1:end-1-21)]);
    yhatALLiw(tt,1+1) = [1,RV1(end-21)]*tempAR.beta;
    
    % 2 HAR
    tempHAR = ols(RV1m(21+2:end),[ones(tt-23-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21)]*tempHAR.beta;

    % 3 HAR-J
    tempHARJ = ols(RV1m(21+2:end),[ones(tt-23-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RJ1(1:end-1-21)]);
    yhatALLiw(tt,1+3) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RJ1(end-21)]*tempHARJ.beta;
    
    % 4 CHAR
    tempCHAR = ols(RV1m(21+2:end),[ones(tt-23-21,1),BPV1d(1:end-1-21),BPV1w(1:end-1-21),BPV1m(1:end-1-21)]);
    yhatALLiw(tt,1+4) = [1,BPV1d(end-21),BPV1w(end-21),BPV1m(end-21)]*tempCHAR.beta;
    
    % 5 SHAR
    tempSHAR = ols(RV1m(21+2:end),[ones(tt-23-21,1),RVp1(1:end-1-21),RVn1(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21)]);
    yhatALLiw(tt,1+5) = [1,RVp1(end-21),RVn1(end-21),RV1w(end-21),RV1m(end-21)]*tempSHAR.beta;
    
    % 6 ARQ
    tempARQ = ols(RV1m(21+2:end),[ones(tt-23-21,1),RV1d(1:end-1-21),RV1d(1:end-1-21).*sqrtRQd(1:end-1-21)]);
    yhatALLiw(tt,1+6) = [1,RV1d(end-21),RV1d(end-21).*sqrtRQd(end-21)]*tempARQ.beta;
    
    % 7 HARQ
    tempHARQ = ols(RV1m(21+2:end),[ones(tt-23-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RV1d(1:end-1-21).*sqrtRQd(1:end-1-21)]);
    yhatALLiw(tt,1+7) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RV1d(end-21).*sqrtRQd(end-21)]*tempHARQ.beta;
    
    % 8 HARQF
    tempHARQF = ols(RV1m(21+2:end),[ones(tt-23-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RV1d(1:end-1-21).*sqrtRQd(1:end-1-21),RV1w(1:end-1-21).*sqrtRQw(1:end-1-21),RV1m(1:end-1-21).*sqrtRQm(1:end-1-21)]);
    yhatALLiw(tt,1+8) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RV1d(end-21).*sqrtRQd(end-21),RV1w(end-21).*sqrtRQw(end-21),RV1m(end-21).*sqrtRQm(end-21)]*tempHARQF.beta;
        
    % 9 HARQF
    tempHARQh = ols(RV1m(21+2:end),[ones(tt-23-21,1),RV1d(1:end-1-21),RV1w(1:end-1-21),RV1m(1:end-1-21),RV1m(1:end-1-21).*sqrtRQm(1:end-1-21)]);
    yhatALLiw(tt,1+9) = [1,RV1d(end-21),RV1w(end-21),RV1m(end-21),RV1m(end-21).*sqrtRQm(end-21)]*tempHARQh.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:9
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countiw(ii) = countiw(ii) + temp;
        end
    end
end

yhatALL = nan(T,1+9,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1001 obs where we have no OOS forecast

table8a = nan(4,9);
table8 = nan(4,9);
for mm=[2,1,(3:9)]
    table8a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table8a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table8a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table8a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window

    table8(:,mm) = table8a(:,mm)./table8a(:,2);
end
clear info;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F','HARQ-h');
info.rnames = char('Table 8','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
info.width = 150;
sprintf('\n\n');
mprint(table8,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 9: OOS (daily) comparisons with other RV measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, the rolling window results
T = length(data1min);
window = 1000;
yhatALLrw = nan(T,1+6);  % will store the realized value and the forecasts
countrw = zeros(6,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = data1min(tt-1000:tt-1,1);
    RQ1 = data1min(tt-1000:tt-1,6);
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    SSRV1d  = data1min(tt-1000:tt-1,2);
    SSRV1w = mean([SSRV1d,mlag(SSRV1d,4,mean(SSRV1d))],2);
    SSRV1m = mean([SSRV1d,mlag(SSRV1d,21,mean(SSRV1d))],2);
    
    TSRV1d = data1min(tt-1000:tt-1,3);
    TSRV1w = mean([TSRV1d,mlag(TSRV1d,4,mean(TSRV1d))],2);
    TSRV1m = mean([TSRV1d,mlag(TSRV1d,21,mean(TSRV1d))],2);
    
    RK1d = data1min(tt-1000:tt-1,4);
    RK1w = mean([RK1d,mlag(RK1d,4,mean(RK1d))],2);
    RK1m = mean([RK1d,mlag(RK1d,21,mean(RK1d))],2);
    
    PARV1d = data1min(tt-1000:tt-1,5);
    PARV1w = mean([PARV1d,mlag(PARV1d,4,mean(PARV1d))],2);
    PARV1m = mean([PARV1d,mlag(PARV1d,21,mean(PARV1d))],2);

    minRV   = min(RV1d);          maxRV = max(RV1d);        meanRV = mean(RV1d);    % will use these values in the "insanity filter" below
    
       RV1 =    RV1(22:end);
      RV1d =   RV1d(22:end);      RV1w =   RV1w(22:end);      RV1m =   RV1m(22:end);
    SSRV1d = SSRV1d(22:end);    SSRV1w = SSRV1w(22:end);    SSRV1m = SSRV1m(22:end);
    TSRV1d = TSRV1d(22:end);    TSRV1w = TSRV1w(22:end);    TSRV1m = TSRV1m(22:end);
      RK1d =   RK1d(22:end);      RK1w =   RK1w(22:end);      RK1m =   RK1m(22:end);
    PARV1d = PARV1d(22:end);    PARV1w = PARV1w(22:end);    PARV1m = PARV1m(22:end);
    
    sqrtRQd = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    
    % getting the realised vol for this day
    yhatALLrw(tt,1) = data1min(tt,1);
    
    % 1 HARQ (on RV)
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;

    % 2 HAR-RV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHARQ.beta;
    
    % 3 HAR-SSRV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),SSRV1d(1:end-1),SSRV1w(1:end-1),SSRV1m(1:end-1)]);
    yhatALLrw(tt,1+3) = [1,SSRV1d(end),SSRV1w(end),SSRV1m(end)]*tempHARQ.beta;

    % 4 HAR-TSRV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),TSRV1d(1:end-1),TSRV1w(1:end-1),TSRV1m(1:end-1)]);
    yhatALLrw(tt,1+4) = [1,TSRV1d(end),TSRV1w(end),TSRV1m(end)]*tempHARQ.beta;

    % 5 HAR-RK
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RK1d(1:end-1),RK1w(1:end-1),RK1m(1:end-1)]);
    yhatALLrw(tt,1+5) = [1,RK1d(end),RK1w(end),RK1m(end)]*tempHARQ.beta;
    
    % 6 HAR-PARV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),PARV1d(1:end-1),PARV1w(1:end-1),PARV1m(1:end-1)]);
    yhatALLrw(tt,1+6) = [1,PARV1d(end),PARV1w(end),PARV1m(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:size(yhatALLrw,2)-1
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+6);  % will store the realized value and the forecasts
countiw = zeros(6,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = data1min(1:tt-1,1);
    RQ1 = data1min(1:tt-1,6);
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    SSRV1d = data1min(1:tt-1,2);
    SSRV1w = mean([SSRV1d,mlag(SSRV1d,4,mean(SSRV1d))],2);
    SSRV1m = mean([SSRV1d,mlag(SSRV1d,21,mean(SSRV1d))],2);
    
    TSRV1d = data1min(1:tt-1,3);
    TSRV1w = mean([TSRV1d,mlag(TSRV1d,4,mean(TSRV1d))],2);
    TSRV1m = mean([TSRV1d,mlag(TSRV1d,21,mean(TSRV1d))],2);
    
    RK1d = data1min(1:tt-1,4);
    RK1w = mean([RK1d,mlag(RK1d,4,mean(RK1d))],2);
    RK1m = mean([RK1d,mlag(RK1d,21,mean(RK1d))],2);
    
    PARV1d = data1min(1:tt-1,5);
    PARV1w = mean([PARV1d,mlag(PARV1d,4,mean(PARV1d))],2);
    PARV1m = mean([PARV1d,mlag(PARV1d,21,mean(PARV1d))],2);

    minRV   = min(data1min(tt-1000:tt-1,1));          maxRV = max(data1min(tt-1000:tt-1,1));        meanRV = mean(data1min(tt-1000:tt-1,1));    % will use these values in the "insanity filter" below
 
       RV1 =    RV1(22:end);
      RV1d =   RV1d(22:end);      RV1w =   RV1w(22:end);      RV1m =   RV1m(22:end);
    SSRV1d = SSRV1d(22:end);    SSRV1w = SSRV1w(22:end);    SSRV1m = SSRV1m(22:end);
    TSRV1d = TSRV1d(22:end);    TSRV1w = TSRV1w(22:end);    TSRV1m = TSRV1m(22:end);
      RK1d =   RK1d(22:end);      RK1w =   RK1w(22:end);      RK1m =   RK1m(22:end);
    PARV1d = PARV1d(22:end);    PARV1w = PARV1w(22:end);    PARV1m = PARV1m(22:end);

    sqrtRQd = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    
    % getting the realised vol for this day
    yhatALLiw(tt,1) = data1min(tt,1);
    
    % 1 HARQ (on RV)
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;

    % 2 HAR-RV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHARQ.beta;
    
    % 3 HAR-SSRV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),SSRV1d(1:end-1),SSRV1w(1:end-1),SSRV1m(1:end-1)]);
    yhatALLiw(tt,1+3) = [1,SSRV1d(end),SSRV1w(end),SSRV1m(end)]*tempHARQ.beta;

    % 4 HAR-TSRV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),TSRV1d(1:end-1),TSRV1w(1:end-1),TSRV1m(1:end-1)]);
    yhatALLiw(tt,1+4) = [1,TSRV1d(end),TSRV1w(end),TSRV1m(end)]*tempHARQ.beta;

    % 5 HAR-RK
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RK1d(1:end-1),RK1w(1:end-1),RK1m(1:end-1)]);
    yhatALLiw(tt,1+5) = [1,RK1d(end),RK1w(end),RK1m(end)]*tempHARQ.beta;
    
    % 6 HAR-PARV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),PARV1d(1:end-1),PARV1w(1:end-1),PARV1m(1:end-1)]);
    yhatALLiw(tt,1+6) = [1,PARV1d(end),PARV1w(end),PARV1m(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:size(yhatALLiw,2)-1
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countiw(ii) = countiw(ii) + temp;
        end
    end
end



yhatALL = nan(T,1+6,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1000 obs where we have no OOS forecast

table9a = nan(4,6);
table9 = nan(4,6);
for mm=1:6
    table9a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table9a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table9a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table9a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window

    table9(:,mm) = table9a(:,mm)./table9a(:,1);
end
clear info;
info.cnames = char('RV','SS-RV','TS-RV','RK','PA-RV');
info.rnames = char('Table 9','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table9(:,2:end),info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 10: HAR vs HARQ using different RV measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, the rolling window results
window = 1000;
yhatALLrw = nan(T,1+10);  % will store the realized value and the forecasts
countrw = zeros(10,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1    = data1min(tt-1000:tt-1,1);
    RQ1    = data1min(tt-1000:tt-1,6);
    RQ1min = data1min(tt-1000:tt-1,6);
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    SSRV1d  = data1min(tt-1000:tt-1,2);
    SSRV1w = mean([SSRV1d,mlag(SSRV1d,4,mean(SSRV1d))],2);
    SSRV1m = mean([SSRV1d,mlag(SSRV1d,21,mean(SSRV1d))],2);
    
    TSRV1d = data1min(tt-1000:tt-1,3);
    TSRV1w = mean([TSRV1d,mlag(TSRV1d,4,mean(TSRV1d))],2);
    TSRV1m = mean([TSRV1d,mlag(TSRV1d,21,mean(TSRV1d))],2);
    
    RK1d = data1min(tt-1000:tt-1,4);
    RK1w = mean([RK1d,mlag(RK1d,4,mean(RK1d))],2);
    RK1m = mean([RK1d,mlag(RK1d,21,mean(RK1d))],2);
    
    PARV1d = data1min(tt-1000:tt-1,5);
    PARV1w = mean([PARV1d,mlag(PARV1d,4,mean(PARV1d))],2);
    PARV1m = mean([PARV1d,mlag(PARV1d,21,mean(PARV1d))],2);

    minRV   = min(RV1d);          maxRV = max(RV1d);        meanRV = mean(RV1d);    % will use these values in the "insanity filter" below
    
       RV1 =    RV1(22:end);
      RV1d =   RV1d(22:end);      RV1w =   RV1w(22:end);      RV1m =   RV1m(22:end);
    SSRV1d = SSRV1d(22:end);    SSRV1w = SSRV1w(22:end);    SSRV1m = SSRV1m(22:end);
    TSRV1d = TSRV1d(22:end);    TSRV1w = TSRV1w(22:end);    TSRV1m = TSRV1m(22:end);
      RK1d =   RK1d(22:end);      RK1w =   RK1w(22:end);      RK1m =   RK1m(22:end);
    PARV1d = PARV1d(22:end);    PARV1w = PARV1w(22:end);    PARV1m = PARV1m(22:end);
    
     sqrtRQd = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    
    % getting the realised vol for this day
    yhatALLrw(tt,1) = data1min(tt,1);
    
    % 1 HAR-RV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLrw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHARQ.beta;

    % 2 HARQ-RV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;

    % 3 HAR-SSRV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),SSRV1d(1:end-1),SSRV1w(1:end-1),SSRV1m(1:end-1)]);
    yhatALLrw(tt,1+3) = [1,SSRV1d(end),SSRV1w(end),SSRV1m(end)]*tempHARQ.beta;

    % 4 HARQ-SSRV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),SSRV1d(1:end-1),SSRV1w(1:end-1),SSRV1m(1:end-1),SSRV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+4) = [1,SSRV1d(end),SSRV1w(end),SSRV1m(end),SSRV1d(end).*sqrtRQd(end)]*tempHARQ.beta;

    % 5 HAR-TSRV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),TSRV1d(1:end-1),TSRV1w(1:end-1),TSRV1m(1:end-1)]);
    yhatALLrw(tt,1+5) = [1,TSRV1d(end),TSRV1w(end),TSRV1m(end)]*tempHARQ.beta;

    % 6 HARQ-TSRV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),TSRV1d(1:end-1),TSRV1w(1:end-1),TSRV1m(1:end-1),TSRV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+6) = [1,TSRV1d(end),TSRV1w(end),TSRV1m(end),TSRV1d(end).*sqrtRQd(end)]*tempHARQ.beta;

    % 7 HAR-RK
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RK1d(1:end-1),RK1w(1:end-1),RK1m(1:end-1)]);
    yhatALLrw(tt,1+7) = [1,RK1d(end),RK1w(end),RK1m(end)]*tempHARQ.beta;
    
    % 8 HARQ-RK
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RK1d(1:end-1),RK1w(1:end-1),RK1m(1:end-1),RK1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+8) = [1,RK1d(end),RK1w(end),RK1m(end),RK1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 9 HAR-PARV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),PARV1d(1:end-1),PARV1w(1:end-1),PARV1m(1:end-1)]);
    yhatALLrw(tt,1+9) = [1,PARV1d(end),PARV1w(end),PARV1m(end)]*tempHARQ.beta;
    
    % 10 HARQ-PARV
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),PARV1d(1:end-1),PARV1w(1:end-1),PARV1m(1:end-1),PARV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+10) = [1,PARV1d(end),PARV1w(end),PARV1m(end),PARV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:size(yhatALLrw,2)-1
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+10);  % will store the realized value and the forecasts
countiw = zeros(10,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = data1min(1:tt-1,1);
    RQ1 = data1min(1:tt-1,6);
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    SSRV1d = data1min(1:tt-1,2);
    SSRV1w = mean([SSRV1d,mlag(SSRV1d,4,mean(SSRV1d))],2);
    SSRV1m = mean([SSRV1d,mlag(SSRV1d,21,mean(SSRV1d))],2);
    
    TSRV1d = data1min(1:tt-1,3);
    TSRV1w = mean([TSRV1d,mlag(TSRV1d,4,mean(TSRV1d))],2);
    TSRV1m = mean([TSRV1d,mlag(TSRV1d,21,mean(TSRV1d))],2);
    
    RK1d = data1min(1:tt-1,4);
    RK1w = mean([RK1d,mlag(RK1d,4,mean(RK1d))],2);
    RK1m = mean([RK1d,mlag(RK1d,21,mean(RK1d))],2);
    
    PARV1d = data1min(1:tt-1,5);
    PARV1w = mean([PARV1d,mlag(PARV1d,4,mean(PARV1d))],2);
    PARV1m = mean([PARV1d,mlag(PARV1d,21,mean(PARV1d))],2);
    
    minRV   = min(data1min(tt-1000:tt-1,1));          maxRV = max(data1min(tt-1000:tt-1,1));        meanRV = mean(data1min(tt-1000:tt-1,1));    % will use these values in the "insanity filter" below
  
    RV1 =    RV1(22:end);
    RV1d =   RV1d(22:end);      RV1w =   RV1w(22:end);      RV1m =   RV1m(22:end);
    SSRV1d = SSRV1d(22:end);    SSRV1w = SSRV1w(22:end);    SSRV1m = SSRV1m(22:end);
    TSRV1d = TSRV1d(22:end);    TSRV1w = TSRV1w(22:end);    TSRV1m = TSRV1m(22:end);
    RK1d =   RK1d(22:end);      RK1w =   RK1w(22:end);      RK1m =   RK1m(22:end);
    PARV1d = PARV1d(22:end);    PARV1w = PARV1w(22:end);    PARV1m = PARV1m(22:end);
    
    sqrtRQd = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    
    % getting the realised vol for this day
    yhatALLiw(tt,1) = data1min(tt,1);
    
    % 1 HAR-RV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLiw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHARQ.beta;
    
    % 2 HARQ-RV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 3 HAR-SSRV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),SSRV1d(1:end-1),SSRV1w(1:end-1),SSRV1m(1:end-1)]);
    yhatALLiw(tt,1+3) = [1,SSRV1d(end),SSRV1w(end),SSRV1m(end)]*tempHARQ.beta;
    
    % 4 HARQ-SSRV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),SSRV1d(1:end-1),SSRV1w(1:end-1),SSRV1m(1:end-1),SSRV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+4) = [1,SSRV1d(end),SSRV1w(end),SSRV1m(end),SSRV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 5 HAR-TSRV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),TSRV1d(1:end-1),TSRV1w(1:end-1),TSRV1m(1:end-1)]);
    yhatALLiw(tt,1+5) = [1,TSRV1d(end),TSRV1w(end),TSRV1m(end)]*tempHARQ.beta;
    
    % 6 HARQ-TSRV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),TSRV1d(1:end-1),TSRV1w(1:end-1),TSRV1m(1:end-1),TSRV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+6) = [1,TSRV1d(end),TSRV1w(end),TSRV1m(end),TSRV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 7 HAR-RK
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RK1d(1:end-1),RK1w(1:end-1),RK1m(1:end-1)]);
    yhatALLiw(tt,1+7) = [1,RK1d(end),RK1w(end),RK1m(end)]*tempHARQ.beta;
    
    % 8 HARQ-RK
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RK1d(1:end-1),RK1w(1:end-1),RK1m(1:end-1),RK1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+8) = [1,RK1d(end),RK1w(end),RK1m(end),RK1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 9 HAR-PARV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),PARV1d(1:end-1),PARV1w(1:end-1),PARV1m(1:end-1)]);
    yhatALLiw(tt,1+9) = [1,PARV1d(end),PARV1w(end),PARV1m(end)]*tempHARQ.beta;
    
    % 10 HARQ-PARV
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),PARV1d(1:end-1),PARV1w(1:end-1),PARV1m(1:end-1),PARV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+10) = [1,PARV1d(end),PARV1w(end),PARV1m(end),PARV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:size(yhatALLiw,2)-1
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countiw(ii) = countiw(ii) + temp;
        end
    end
end


yhatALL = nan(T,1+10,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1000 obs where we have no OOS forecast

table10a = nan(4,10);
table10 = nan(4,10);
for mm=1:10
    table10a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table10a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table10a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table10a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window
end

clear info;
info.cnames = char('RV','SS-RV','TS-RV','RK','PA-RV');
info.rnames = char('Table 10','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table10a(:,2:2:end)./table10a(:,1:2:end),info)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 11: Different IQ estimators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tables 11 and 12 use the SP500 index data again.

% first, the rolling window results
T = length(RV);
window = 1000;
yhatALLrw = nan(T,1+6);  % will store the realized value and the forecasts
countrw = zeros(6,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(tt-1000:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(tt-1000:tt-1);
    
    minRV = min(RV1);  % will use these values in the "insanity filter" below
    maxRV = max(RV1);
    meanRV = mean(RV1);
 
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    RV1 = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    TPQ1     = sqrt(TPQ(tt-1000+21:tt-1));  
    MedRQ1   = sqrt(MedRQ(tt-1000+21:tt-1));
    TrRQ1   = sqrt(TrRQ(tt-1000+21:tt-1));
    RQ15min1 = sqrt(RQ15min(tt-1000+21:tt-1));
    RQboot1  = sqrt(RQboot(tt-1000+21:tt-1));
    
    % getting the realised vol for this day
    yhatALLrw(tt,1) = RV(tt);
    
    % 1 RQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 2 TPQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*TPQ1(1:end-1)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*TPQ1(end)]*tempHARQ.beta;

    % 3 MedRQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*MedRQ1(1:end-1)]);
    yhatALLrw(tt,1+3) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*MedRQ1(end)]*tempHARQ.beta;

    % 4 TrRQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*TrRQ1(1:end-1)]);
    yhatALLrw(tt,1+4) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*TrRQ1(end)]*tempHARQ.beta;
    
    % 5 RQ15min
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*RQ15min1(1:end-1)]);
    yhatALLrw(tt,1+5) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*RQ15min1(end)]*tempHARQ.beta;
    
    % 6 RQboot
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*RQboot1(1:end-1)]);
    yhatALLrw(tt,1+6) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*RQboot1(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:6
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+6);  % will store the realized value and the forecasts
countiw = zeros(6,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(1:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(1:tt-1);
    
    minRV = min(RV(tt-1000:tt-1));  % will use these values in the "insanity filter" below
    maxRV = max(RV(tt-1000:tt-1));
    meanRV = mean(RV(tt-1000:tt-1));
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    RV1  = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd  = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    TPQ1     = sqrt(TPQ(1+21:tt-1));  
    MedRQ1   = sqrt(MedRQ(1+21:tt-1));
    TrRQ1    = sqrt(TrRQ(1+21:tt-1));
    RQ15min1 = sqrt(RQ15min(1+21:tt-1));
    RQboot1  = sqrt(RQboot(1+21:tt-1));
    
    % getting the realised vol for this day
    yhatALLiw(tt,1) = RV(tt);
    
   % 1 RQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 2 TPQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*TPQ1(1:end-1)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*TPQ1(end)]*tempHARQ.beta;

    % 3 MedRQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*MedRQ1(1:end-1)]);
    yhatALLiw(tt,1+3) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*MedRQ1(end)]*tempHARQ.beta;

    % 4 TrRQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*TrRQ1(1:end-1)]);
    yhatALLiw(tt,1+4) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*TrRQ1(end)]*tempHARQ.beta;
    
    % 5 RQ15min
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*RQ15min1(1:end-1)]);
    yhatALLiw(tt,1+5) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*RQ15min1(end)]*tempHARQ.beta;
    
    % 6 RQboot
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*RQboot1(1:end-1)]);
    yhatALLiw(tt,1+6) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*RQboot1(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:6
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countiw(ii) = countiw(ii) + temp;
        end
    end
end


yhatALL = nan(T,1+6,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1000 obs where we have no OOS forecast

table11a = nan(4,6);
table11 = nan(4,6);
for mm=1:6
    table11a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table11a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table11a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table11a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window

    table11(:,mm) = table11a(:,mm)./table11a(:,1);
end
clear info;
info.cnames = char('RQ','TPQ','MedRQ','TrRQ','RQ15min','Boot');
info.rnames = char('Table 11','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table11,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 12: Different HARQ specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% first, the rolling window results
window = 1000;
yhatALLrw = nan(T,1+7);  % will store the realized value and the forecasts
countrw = zeros(7,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(tt-1000:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(tt-1000:tt-1);
    
    minRV = min(RV1);  % will use these values in the "insanity filter" below
    maxRV = max(RV1);
    meanRV = mean(RV1);
 
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    RV1  = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd    = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    invsqrtRQd = 1./sqrt(RQ1(22:end));  
    RQd        = RQ1(22:end);
    invRQd     = 1./RQ1(22:end);
    logRQd     = log(RQ1(22:end));
    
    % getting the realised vol for this day
    yhatALLrw(tt,1) = RV(tt);
    
    % 1 RQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*RQd(1:end-1)]);
    yhatALLrw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*RQd(end)]*tempHARQ.beta;
    
    % 2 sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;

    % 3 1/sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*invsqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+3) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*invsqrtRQd(end)]*tempHARQ.beta;
    
    % 4 1/RQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*invRQd(1:end-1)]);
    yhatALLrw(tt,1+4) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*invRQd(end)]*tempHARQ.beta;
    
    % 5 log(RQ)
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*logRQd(1:end-1)]);
    yhatALLrw(tt,1+5) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*logRQd(end)]*tempHARQ.beta;
    
    % 6 HAR with level of sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+6) = [1,RV1d(end),RV1w(end),RV1m(end),sqrtRQd(end)]*tempHARQ.beta;
    
    % 7 HARQ with level of sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1),sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+7) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end),sqrtRQd(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:7
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+7);  % will store the realized value and the forecasts
countiw = zeros(7,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(1:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(1:tt-1);
    
    minRV = min(RV(tt-1000:tt-1));  % will use these values in the "insanity filter" below
    maxRV = max(RV(tt-1000:tt-1));
    meanRV = mean(RV(tt-1000:tt-1));
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    RV1 = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);

    sqrtRQd = sqrt(RQ1(22:end));  % won't bother with de-meaning here, as fit is unaffected
    invsqrtRQd = 1./sqrt(RQ1(22:end));
    RQd = RQ1(22:end);
    invRQd = 1./RQ1(22:end);
    logRQd = log(RQ1(22:end));
    
    % getting the realised vol for this day
    yhatALLiw(tt,1) = RV(tt);
    
    % 1 RQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*RQd(1:end-1)]);
    yhatALLiw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*RQd(end)]*tempHARQ.beta;
    
    % 2 sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;

    % 3 1/sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*invsqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+3) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*invsqrtRQd(end)]*tempHARQ.beta;
    
    % 4 1/RQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*invRQd(1:end-1)]);
    yhatALLiw(tt,1+4) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*invRQd(end)]*tempHARQ.beta;
    
    % 5 log(RQ)
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*logRQd(1:end-1)]);
    yhatALLiw(tt,1+5) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*logRQd(end)]*tempHARQ.beta;
    
    % 6 HAR with level of sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+6) = [1,RV1d(end),RV1w(end),RV1m(end),sqrtRQd(end)]*tempHARQ.beta;
    
    % 7 HARQ with level of sqrt(RQ)
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1),sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+7) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end),sqrtRQd(end)]*tempHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:7
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countiw(ii) = countiw(ii) + temp;
        end
    end
end



yhatALL = nan(T,1+7,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1000 obs where we have no OOS forecast

table12a = nan(4,7);
table12 = nan(4,7);
for mm=[2,1,(3:7)]
    table12a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table12a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table12a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table12a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window

    table12(:,mm) = table12a(:,mm)./table12a(:,2);
end
clear info;
info.cnames = char('RQ','RQ^1/2','RQ^-1/2','RQ^-1','Log RQ','HAR','HARQ');
info.rnames = char('Table 12','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table12,info)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 13: Alternative Q-model in sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RVd = RV;
RVw = mean([RV,mlag(RV,4,mean(RV))],2);
RVm = mean([RV,mlag(RV,21,mean(RV))],2);  % note these are not lagged, so need to lag one period if using in a forecasting regression

BPVd = BPV;
BPVw = mean([BPV,mlag(BPV,4,mean(BPV))],2);
BPVm = mean([BPV,mlag(BPV,21,mean(BPV))],2);  

sqrtRQdemeanD = sqrt(RQ) - mean(sqrt(RQ));  % we use this in our "HARQ" model
sqrtTPQdemeanD = sqrt(TPQ)-mean(sqrt(TPQ));
% HAR-J
tempHARJ = hwhite(RV(23:end),[ones(T-22,1),RVd(22:end-1),RVw(22:end-1), RVm(22:end-1),RJ(22:end-1)]);
hhatHARJ = [ones(T-22,1),RVd(22:end-1),RVw(22:end-1), RVm(22:end-1),RJ(22:end-1)]*tempHARJ.beta;  % dropping first 21 obs to make comparable with HAR

% HARQ-J
tempHARQJ = hwhite(RV(23:end),[ones(T-22,1),RVd(22:end-1),RVw(22:end-1), RVm(22:end-1),RJ(22:end-1),RVd(22:end-1).*sqrtRQdemeanD(22:end-1)]);
hhatHARQJ = [ones(T-22,1),RVd(22:end-1),RVw(22:end-1), RVm(22:end-1),RJ(22:end-1),RVd(22:end-1).*RQ(22:end-1)]*tempHARQJ.beta;

% CHAR
tempCHAR = hwhite(RV(23:end),[ones(T-22,1),BPVd(22:end-1),BPVw(22:end-1),BPVm(22:end-1)]);
hhatCHAR = [ones(T-22,1),BPVd(22:end-1),BPVw(22:end-1),BPVm(22:end-1)]*tempCHAR.beta;

% CHARQ
tempCHARQ = hwhite(RV(23:end),[ones(T-22,1),BPVd(22:end-1),BPVw(22:end-1),BPVm(22:end-1),BPVd(22:end-1).*sqrtTPQdemeanD(22:end-1)]);
hhatCHARQ = [ones(T-22,1),BPVd(22:end-1),BPVw(22:end-1),BPVm(22:end-1),BPVd(22:end-1).*sqrtTPQdemeanD(22:end-1)]*tempCHARQ.beta;

% SHAR
tempSHAR = hwhite(RV(23:end),[ones(T-22,1),RVw(22:end-1),RVm(22:end-1),RVp(22:end-1),RVn(22:end-1)]);
hhatSHAR = [ones(T-22,1),RVw(22:end-1),RVm(22:end-1),RVp(22:end-1),RVn(22:end-1)]*tempSHAR.beta;

% SHARQ
tempSHARQ = hwhite(RV(23:end),[ones(T-22,1),RVw(22:end-1),RVm(22:end-1),RVp(22:end-1),RVn(22:end-1),RVp(22:end-1).*sqrtRQdemeanD(22:end-1), RVd(22:end-1).*sqrtRQdemeanD(22:end-1)]);
hhatSHARQ = [ones(T-22,1),RVw(22:end-1),RVm(22:end-1),RVp(22:end-1),RVn(22:end-1),RVp(22:end-1).*sqrtRQdemeanD(22:end-1), RVd(22:end-1).*sqrtRQdemeanD(22:end-1)]*tempSHARQ.beta;

table13 = nan(10*2+3,6);
table13([1,3,5,7,9],1) = tempHARJ.beta;
table13([2,4,6,8,10],1) = tempHARJ.se;
table13([1,3,5,7,9,15],2) = tempHARQJ.beta;
table13([2,4,6,8,10,16],2) = tempHARQJ.se;

table13([1,3,5,7],3) = tempCHAR.beta;
table13([2,4,6,8],3) = tempCHAR.se;
table13([1,3,5,7,15],4) = tempCHARQ.beta;
table13([2,4,6,8,16],4) = tempCHARQ.se;

table13([1,5,7,11,13],5) = tempSHAR.beta;
table13([2,6,8,12,14],5) = tempSHAR.se;
table13([1,5,7,11,13,17,19],6) = tempSHARQ.beta;
table13([2,6,8,12,14,18,20],6) = tempSHARQ.se;

hhatALLis = [hhatHARJ,hhatHARQJ,hhatCHAR,hhatCHARQ,hhatSHAR,hhatSHARQ];
for ii=1:size(hhatALLis,2);
    table13(21,ii) = 1-cov(RV(23:end)-hhatALLis(:,ii))/cov(RV(23:end));            % r2
    table13(22,ii) = mean( (RV(23:end)-hhatALLis(:,ii)).^2);                              % MSE
    table13(23,ii) = mean( RV(23:end)./hhatALLis(:,ii) - log(RV(23:end)./hhatALLis(:,ii)) -1 );     % QLIKE
end

clear info;
info.cnames = char('HAR-J', 'HARQ-J', 'CHAR', 'CHARQ', 'SHAR', 'SHARQ');
info.rnames = char('Table 13','b0','s.e.','b1','s.e.','b2','s.e.','b3','s.e.','bJ','s.e.','b1p','s.e.','b1n','s.e.','b1Q','s.e.','b1Qp','s.e.','b1Qn','s.e.','R2','MSE','QLIKE');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table13,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 14: out-of-sample alternative models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first, the rolling window results
window = 1000;
yhatALLrw = nan(T,1+8);  % will store the realized value and the forecasts: AR, HAR, HAR-J, CHAR, SHAR, ARQ, HARQ, HARQF
countrw = zeros(8,1);


for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(tt-1000:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(tt-1000:tt-1);
    BPV1 = BPV(tt-1000:tt-1);
    TPQ1 = TPQ(tt-1000:tt-1);
    
    minRV = min(RV1);  % will use these values in the "insanity filter" below
    maxRV = max(RV1);
    meanRV = mean(RV1);
 
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1,4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    
    RV1 = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);
    
    BPV1 = BPV1(22:end);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtTPQd = sqrt(TPQ1);
  
    sqrtRQd = sqrtRQd(22:end);
    sqrtTPQd = sqrtTPQd(22:end);
    
    RJ1  = RJ(tt-1000+21:tt-1);
    RVp1 = RVp(tt-1000+21:tt-1);
    RVn1 = RVn(tt-1000+21:tt-1);
    
    % getting the realised vol for this day
    yhatALLrw(tt,1) = RV(tt);
    
    % 1 HAR
    tempHAR = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLrw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHAR.beta;
    
    % 2 HARQ
    tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 3 HAR-J
    tempHARJ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1), RV1m(1:end-1),RJ1(1:end-1)]);
    yhatALLrw(tt,1+3) = [1,RV1d(end),RV1w(end), RV1m(end),RJ1(end)]*tempHARJ.beta;  

    % 4 HARQ-J
    tempHARQJ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1), RV1m(1:end-1),RJ1(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+4) = [1,RV1d(end),RV1w(end), RV1m(end),RJ1(end),RV1d(end).*sqrtRQd(end)]*tempHARQJ.beta;

    % 5 CHAR
    tempCHAR = ols(RV1(2:end),[ones(window-22,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1)]);
    yhatALLrw(tt,1+5) = [1,BPV1d(end),BPV1w(end),BPV1m(end)]*tempCHAR.beta;

    % 6 CHARQ
    tempCHARQ = ols(RV1(2:end),[ones(window-22,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1),BPV1d(1:end-1).*sqrtTPQd(1:end-1)]);
    yhatALLrw(tt,1+6) = [1,BPV1d(end),BPV1w(end),BPV1m(end),BPV1d(end).*sqrtTPQd(end)]*tempCHARQ.beta;

    % 7 SHAR
    tempSHAR = ols(RV1(2:end),[ones(window-22,1),RV1w(1:end-1),RV1m(1:end-1),RVp1(1:end-1),RVn1(1:end-1)]);
    yhatALLrw(tt,1+7) = [1,RV1w(end),RV1m(end),RVp1(end),RVn1(end)]*tempSHAR.beta;

    % 8 SHARQ
    tempSHARQ = ols(RV1(2:end),[ones(window-22,1),RVp1(1:end-1),RVn1(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RVp1(1:end-1).*sqrtRQd(1:end-1), RVn1(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLrw(tt,1+8) = [1,RV1w(end),RV1m(end),RVp1(end),RVn1(end),RVp1(end).*sqrtRQd(end), RVn1(end).*sqrtRQd(end)]*tempSHARQ.beta;
    
    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:8
            [yhatALLrw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


% second, the expanding window results
yhatALLiw = nan(T,1+8);  % will store the realized value and the forecasts
countiw = zeros(8,1);

for tt=window+1:T;
    % pulling out the data and getting it ready for the models
    RV1 = RV(1:tt-1); % this is the data available to forecast RV(tt)
    RQ1 = RQ(1:tt-1);
    TPQ1 = TPQ(1:tt-1);
    
    minRV = min(RV(tt-1000:tt-1));  % will use these values in the "insanity filter" below
    maxRV = max(RV(tt-1000:tt-1));
    meanRV = mean(RV(tt-1000:tt-1));
    
    RV1d = RV1;
    RV1w = mean([RV1,mlag(RV1,4,mean(RV1))],2);
    RV1m = mean([RV1,mlag(RV1,21,mean(RV1))],2);
    
    RV1  = RV1(22:end);
    RV1d = RV1d(22:end);
    RV1w = RV1w(22:end);
    RV1m = RV1m(22:end);
    
    RJ1  = RJ(1+21:tt-1);
    RVp1 = RVp(1+21:tt-1);
    RVn1 = RVn(1+21:tt-1);

    BPV1  = BPV(1:tt-1);
    BPV1d = BPV1;
    BPV1w = mean([BPV1,mlag(BPV1,4,mean(BPV1))],2);
    BPV1m = mean([BPV1,mlag(BPV1,21,mean(BPV1))],2);
    
    BPV1 = BPV1(22:end);
    BPV1d = BPV1d(22:end);
    BPV1w = BPV1w(22:end);
    BPV1m = BPV1m(22:end);

    sqrtRQd = sqrt(RQ1);  % won't bother with de-meaning here, as fit is unaffected
    sqrtTPQd = sqrt(TPQ1);
    
    sqrtRQd = sqrtRQd(22:end);
    sqrtTPQd = sqrtTPQd(22:end);
  
    % getting the realised vol for this day
    yhatALLiw(tt,1) = RV(tt);
    
    % 1 HAR
    tempHAR = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
    yhatALLiw(tt,1+1) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHAR.beta;
    
    % 2 HARQ
    tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+2) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
    
    % 3 HAR-J
    tempHARJ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1), RV1m(1:end-1),RJ1(1:end-1)]);
    yhatALLiw(tt,1+3) = [1,RV1d(end),RV1w(end), RV1m(end),RJ1(end)]*tempHARJ.beta;  

    % 4 HARQ-J
    tempHARQJ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1), RV1m(1:end-1),RJ1(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+4) = [1,RV1d(end),RV1w(end), RV1m(end),RJ1(end),RV1d(end).*sqrtRQd(end)]*tempHARQJ.beta;

    % 5 CHAR
    tempCHAR = ols(RV1(2:end),[ones(tt-23,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1)]);
    yhatALLiw(tt,1+5) = [1,BPV1d(end),BPV1w(end),BPV1m(end)]*tempCHAR.beta;

    % 6 CHARQ
    tempCHARQ = ols(RV1(2:end),[ones(tt-23,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1),BPV1d(1:end-1).*sqrtTPQd(1:end-1)]);
    yhatALLiw(tt,1+6) = [1,BPV1d(end),BPV1w(end),BPV1m(end),BPV1d(end).*sqrtTPQd(end)]*tempCHARQ.beta;

    % 7 SHAR
    tempSHAR = ols(RV1(2:end),[ones(tt-23,1),RV1w(1:end-1),RV1m(1:end-1),RVp1(1:end-1),RVn1(1:end-1)]);
    yhatALLiw(tt,1+7) = [1,RV1w(end),RV1m(end),RVp1(end),RVn1(end)]*tempSHAR.beta;

    % 8 SHARQ
    tempSHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1w(1:end-1),RV1m(1:end-1),RVp1(1:end-1),RVn1(1:end-1),RVp1(1:end-1).*sqrtRQd(1:end-1), RVn1(1:end-1).*sqrtRQd(1:end-1)]);
    yhatALLiw(tt,1+8) = [1,RV1w(end),RV1m(end),RVp1(end),RVn1(end),RVp1(end).*sqrtRQd(end), RVn1(end).*sqrtRQd(end)]*tempSHARQ.beta;

    % applying an insanity filter, if requested
    if insanityfilter==1
        for ii=1:8
            [yhatALLiw(tt,1+ii),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii),minRV,maxRV,meanRV);
            countrw(ii) = countrw(ii) + temp;
        end
    end
end


yhatALL = nan(T,1+8,2);
yhatALL(:,:,1) = yhatALLrw;
yhatALL(:,:,2) = yhatALLiw;
yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1000 obs where we have no OOS forecast

table14a = nan(4,8);
table14 = nan(4,4);
for mm=1:8
    table14a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
    table14a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
    table14a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
    table14a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window
end
table14(:,1) = table14a(:,2)./table14a(:,1);
table14(:,2) = table14a(:,4)./table14a(:,3);
table14(:,3) = table14a(:,6)./table14a(:,5);
table14(:,4) = table14a(:,8)./table14a(:,7);
clear info;
info.cnames = char('HARQ', 'HARQ-J', 'CHARQ', 'SHARQ');
info.rnames = char('Table 14','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
mprint(table14,info)

toc