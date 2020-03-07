% Replicating individual stock results in:
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


% This code takes about 95 mins to run. Almost all of this is
% in the rolling and expanding window estimations, each of which takes
% about 12-14 minutes to run (across the 27 stocks).



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% locating the data file names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

%loading in the raw data files
data_dir = 'E:\BPQ2016_Data';   % change to location where individual data is stored

dt     = dir([data_dir,'\Data_5min\*_RV_5min.xlsx']);  % getting paths to all 5min data files. % change to location where individual data is stored
dt1min = dir([data_dir,'\Data_1min\*_RV_1min.xlsx']);  % getting paths to all 1min data files. % change to location where individual data is stored
num_stocks = length(dt);  % number of files to go through:  27

format short g

insanityfilter = 1;  % =1 if want filter, equal 0 if not

stock_str = {'AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 2: summary statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table2 = nan(num_stocks,6);
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);
    RV  = data(:,1);
    RQ  = data(:,2);
    T = length(RV);
    
    table2(ss,1:4) = [min(RV),mean(RV),median(RV),max(RV)];
    table2(ss,5) = corrcoef12(RV(2:end),RV(1:end-1));
    % ARQ
    sqrtRQdemean = sqrt(RQ) - mean(sqrt(RQ));  % we use this in our "HARQ" model
    tempARQ = hwhite(RV(23:end),[ones(T-22,1),RV(22:end-1),RV(22:end-1).*sqrtRQdemean(22:end-1)]);
    table2(ss,6) = tempARQ.beta(2);
end

clear info;
info.rnames = char('Table 2', 'AXP', 'BA', 'CAT', 'CSCO', 'CVX', 'DD', 'DIS', 'GE', 'HD','IBM', 'INTC','JNJ', 'JPM','KO', 'MCD', 'MMM', 'MRK', 'MSFT', 'NKE', 'PFE', 'PG','TRV', 'UNH', 'UTX', 'VZ', 'WMT', 'XOM');
info.cnames = char('Min','Mean','Median','Max','AR','ARQ');
info.fmt    = '%10.3f';
sprintf('\n\n');
mprint(table2,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 3: in-sample results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table3 = zeros(3,5);

QLIKE = zeros(num_stocks,5);
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);
    RV  = data(:,1);
    RQ  = data(:,2);
    T = length(RV);
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
    
    hhatALLis = [hhatAR,hhatHAR,hhatARQ,hhatHARQ,hhatHARQF];
    for ii=1:size(hhatALLis,2)
        table3(1,ii) = table3(1,ii) + (1 - cov(RV(23:end)-hhatALLis(:,ii))/cov(RV(23:end))) /num_stocks;            % r2
        table3(2,ii) = table3(2,ii) + (mean( (RV(23:end)-hhatALLis(:,ii)).^2 ))/num_stocks;                              % MSE
        table3(3,ii) = table3(3,ii) + (mean( (RV(23:end)./hhatALLis(:,ii) - log(RV(23:end)./hhatALLis(:,ii)) -1 ).*(hhatALLis(:,ii)>=0) ) )/num_stocks;     % QLIKE. Not pretty, but ignore negative fitted values
    end
end

clear info;
info.cnames = char('AR','HAR','ARQ','HARQ','HARQF');
info.rnames = char('Table 3','R2','MSE','QLIKE');
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
countRWIW = nan(8,num_stocks,2);

yhatALLrw = nan(T,1+8, num_stocks);
yhatALLiw = nan(T,1+8, num_stocks);

table2s = nan(4,8,num_stocks);
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);  % loading the data for the ss^th stock
    % data is:  RV	BPV	MedRV	TrRV	RV_min	RV_plus	RQ	TPQ	MedRQ	TrRQ	RQ_min	RQ_plus	J	RQ_slow	TPQ_slow	MedRQ_slow	TrRQ_slow	RQ_min_slow	RQ_plus_slow	RQ_boot_1	RQ_boot_2
    RV  = data(:,1);
    RQ  = data(:,2);
    RJ  = data(:,3);
    BPV = data(:,4);
    RVn = data(:,5);
    RVp = data(:,6);
    
    T = length(RV);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % out-of-sample results: replicating part of Table 4 in the paper
    
    % first, the rolling window results
    window = 1000;
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
        yhatALLrw(tt,1,ss) = RV(tt);
        
        % 1 AR
        tempAR = ols(RV1(2:end),[ones(window-22,1),RV1(1:end-1)]);
        yhatALLrw(tt,1+1,ss) = [1,RV1(end)]*tempAR.beta;
        
        % 2 HAR
        tempHAR = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
        yhatALLrw(tt,1+2,ss) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHAR.beta;
        
        % 3 HAR-J
        tempHARJ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RJ1(1:end-1)]);
        yhatALLrw(tt,1+3,ss) = [1,RV1d(end),RV1w(end),RV1m(end),RJ1(end)]*tempHARJ.beta;
        
        % 4 CHAR
        tempCHAR = ols(RV1(2:end),[ones(window-22,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1)]);
        yhatALLrw(tt,1+4,ss) = [1,BPV1d(end),BPV1w(end),BPV1m(end)]*tempCHAR.beta;
        
        % 5 SHAR
        tempSHAR = ols(RV1(2:end),[ones(window-22,1),RVp1(1:end-1),RVn1(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
        yhatALLrw(tt,1+5,ss) = [1,RVp1(end),RVn1(end),RV1w(end),RV1m(end)]*tempSHAR.beta;
        
        % 6 ARQ
        tempARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
        yhatALLrw(tt,1+6,ss) = [1,RV1d(end),RV1d(end).*sqrtRQd(end)]*tempARQ.beta;
        
        % 7 HARQ
        tempHARQ = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
        yhatALLrw(tt,1+7,ss) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
        
        % 8 HARQF
        tempHARQF = ols(RV1(2:end),[ones(window-22,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1),RV1w(1:end-1).*sqrtRQw(1:end-1),RV1m(1:end-1).*sqrtRQm(1:end-1)]);
        yhatALLrw(tt,1+8,ss) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end),RV1w(end).*sqrtRQw(end),RV1m(end).*sqrtRQm(end)]*tempHARQF.beta;
        
        % applying an insanity filter, if requested
        if insanityfilter==1
            for ii=1:8
                [yhatALLrw(tt,1+ii,ss),temp] = volatility_insanity_filter(yhatALLrw(tt,1+ii,ss),minRV,maxRV,meanRV);
                countrw(ii) = countrw(ii) + temp;
            end
        end
    end
    countRWIW(:,ss,1) = countrw;
    
    % second, the expanding window results
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
        yhatALLiw(tt,1,ss) = RV(tt);
        
        % 1 AR
        tempAR = ols(RV1(2:end),[ones(tt-23,1),RV1(1:end-1)]);
        yhatALLiw(tt,1+1,ss) = [1,RV1(end)]*tempAR.beta;
        
        % 2 HAR
        tempHAR = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
        yhatALLiw(tt,1+2,ss) = [1,RV1d(end),RV1w(end),RV1m(end)]*tempHAR.beta;
        
        % 3 HAR-J
        tempHARJ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RJ1(1:end-1)]);
        yhatALLiw(tt,1+3,ss) = [1,RV1d(end),RV1w(end),RV1m(end),RJ1(end)]*tempHARJ.beta;
        
        % 4 CHAR
        tempCHAR = ols(RV1(2:end),[ones(tt-23,1),BPV1d(1:end-1),BPV1w(1:end-1),BPV1m(1:end-1)]);
        yhatALLiw(tt,1+4,ss) = [1,BPV1d(end),BPV1w(end),BPV1m(end)]*tempCHAR.beta;
        
        % 5 SHAR
        tempSHAR = ols(RV1(2:end),[ones(tt-23,1),RVp1(1:end-1),RVn1(1:end-1),RV1w(1:end-1),RV1m(1:end-1)]);
        yhatALLiw(tt,1+5,ss) = [1,RVp1(end),RVn1(end),RV1w(end),RV1m(end)]*tempSHAR.beta;
        
        % 6 ARQ
        tempARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
        yhatALLiw(tt,1+6,ss) = [1,RV1d(end),RV1d(end).*sqrtRQd(end)]*tempARQ.beta;
        
        % 7 HARQ
        tempHARQ = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1)]);
        yhatALLiw(tt,1+7,ss) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end)]*tempHARQ.beta;
        
        % 8 HARQF
        tempHARQF = ols(RV1(2:end),[ones(tt-23,1),RV1d(1:end-1),RV1w(1:end-1),RV1m(1:end-1),RV1d(1:end-1).*sqrtRQd(1:end-1),RV1w(1:end-1).*sqrtRQw(1:end-1),RV1m(1:end-1).*sqrtRQm(1:end-1)]);
        yhatALLiw(tt,1+8,ss) = [1,RV1d(end),RV1w(end),RV1m(end),RV1d(end).*sqrtRQd(end),RV1w(end).*sqrtRQw(end),RV1m(end).*sqrtRQm(end)]*tempHARQF.beta;
        
        % applying an insanity filter, if requested
        if insanityfilter==1
            for ii=1:8
                [yhatALLiw(tt,1+ii,ss),temp] = volatility_insanity_filter(yhatALLiw(tt,1+ii,ss),minRV,maxRV,meanRV);
                countiw(ii) = countiw(ii) + temp;
            end
        end
    end
    countRWIW(:,ss,2) = countiw;
    
    yhatALL = nan(T,1+8,2);
    yhatALL(:,:,1) = yhatALLrw(:,:,ss);
    yhatALL(:,:,2) = yhatALLiw(:,:,ss);
    yhatALL = yhatALL(window+1:end,:,:); % dropping the first 1002 obs where we have no OOS forecast
    
    table2 = nan(4,8);
    for mm=1:8
        table2(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
        table2(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
        
        table2(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
        table2(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window
    end
    table2s(:,:,ss) = table2;
    
end  % takes about 12 mins for all 27 stocks

yhatALLiw = yhatALLiw(window+1:end,:,:); % dropping the first 1002 obs where we have no OOS forecast
yhatALLrw = yhatALLrw(window+1:end,:,:); % dropping the first 1002 obs where we have no OOS forecast
table2m = mean(table2s,3);

clear info;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F');
info.rnames = char('*Mean avg loss*','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
info.width = 120;
%mprint(table2m,info)

table3m = nan(4,8,2);
for ii=1:8;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table2s(:,ii,ss)./table2s(:,2,ss);
    end
    table3m(:,ii,1) =   mean(temp,2);
    table3m(:,ii,2) = median(temp,2);
end
table33 = nan(8,8);
table33(1:2:end,:) = table3m(:,:,1);
table33(2:2:end,:) = table3m(:,:,2);

info.rnames = char('Table 4','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table33,info)

if 0  % the code below prints out more details on the individual stock results
    clear info;
    info.width = 120;
    info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F');
    info.rnames = char('*MSE-RW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table2s(1,:,:))',info)
    info.rnames = char('*MSE-IW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table2s(2,:,:))',info)
    info.rnames = char('*QLIKE-RW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table2s(3,:,:))',info)
    info.rnames = char('*QLIKE-IW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table2s(4,:,:))',info)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 5: stratified OOS losses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table5 = nan(8+8,8);
table5a = nan(4+4,8,ss);
table5b = nan(4+4,8,ss);

for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);  % loading the data for the ss^th stock
    % data is:  RV	BPV	MedRV	TrRV	RV_min	RV_plus	RQ	TPQ	MedRQ	TrRQ	RQ_min	RQ_plus	J	RQ_slow	TPQ_slow	MedRQ_slow	TrRQ_slow	RQ_min_slow	RQ_plus_slow	RQ_boot_1	RQ_boot_2
    RV  = data(:,1);
    RQ  = data(:,2);
    RJ  = data(:,3);
    BPV = data(:,4);
    RVn = data(:,5);
    RVp = data(:,6);
    
    T = length(RV);
    
    RQ1 = RQ(window:end-1); % pulling out RQ for the OOS period.
    
    for mm=[2,1,(3:8)]
        temp124 = find(RQ1<=quantile(RQ1,0.95)); % days with low RQ
        temp125 = find(RQ1>quantile(RQ1,0.95));  % days with high RQ
        
        table5a(1,mm,ss) = mean( (yhatALLrw(temp124,1,ss)-yhatALLrw(temp124,mm+1,ss)).^2 );% MSE - rolling window
        table5a(2,mm,ss) = mean( (yhatALLiw(temp124,1,ss)-yhatALLiw(temp124,mm+1,ss)).^2 );  % MSE - expanding window
        table5a(3,mm,ss) = mean( yhatALLrw(temp124,1,ss)./yhatALLrw(temp124,mm+1,ss) - log(yhatALLrw(temp124,1,ss)./yhatALLrw(temp124,mm+1,ss)) -1 );  % QLIKE - rolling window
        table5a(4,mm,ss) = mean( yhatALLiw(temp124,1,ss)./yhatALLiw(temp124,mm+1,ss) - log(yhatALLiw(temp124,1,ss)./yhatALLiw(temp124,mm+1,ss)) -1 );  % QLIKE - expanding window
        
        table5a(5,mm,ss) = mean( (yhatALLrw(temp125,1,ss)-yhatALLrw(temp125,mm+1,ss)).^2 );  % MSE - rolling window
        table5a(6,mm,ss) = mean( (yhatALLiw(temp125,1,ss)-yhatALLiw(temp125,mm+1,ss)).^2 );  % MSE - expanding window
        table5a(7,mm,ss) = mean( yhatALLrw(temp125,1,ss)./yhatALLrw(temp125,mm+1,ss) - log(yhatALLrw(temp125,1,ss)./yhatALLrw(temp125,mm+1,ss)) -1 );  % QLIKE - rolling window
        table5a(8,mm,ss) = mean( yhatALLiw(temp125,1,ss)./yhatALLiw(temp125,mm+1,ss) - log(yhatALLiw(temp125,1,ss)./yhatALLiw(temp125,mm+1,ss)) -1 );  % QLIKE - expanding window
        
        table5b(:,mm,ss) = table5a(:,mm,ss)./table5a(:,2,ss);
    end
end

for ii=1:8
    table5((2*ii-1),:)    =   mean(table5b(ii,:,:),3);
    table5(2*ii,:)        =   median(table5b(ii,:,:),3);
end
clear info;
info.width = 120;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F');
info.rnames = char('.','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
info.fmt    = '%10.4f';
sprintf('\n\nTable 5: Bottom 95% RQ')
mprint(table5(1:8,:),info)
sprintf('Table 5: Top 5% RQ')
mprint(table5(9:16,:),info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 6: in-sample weekly and monthly parameter estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Table 6 only used data for the S&P 500, no individual stock data.
sprintf('\n\nTable 6: In-sample weekly and monthly parameter estimates. \nOnly used data for the S&P 500, no individual stock data.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 7: OOS weekly estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table7s = nan(4,9,num_stocks);
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);  % loading the data for the ss^th stock
    % data is:  RV	BPV	MedRV	TrRV	RV_min	RV_plus	RQ	TPQ	MedRQ	TrRQ	RQ_min	RQ_plus	J	RQ_slow	TPQ_slow	MedRQ_slow	TrRQ_slow	RQ_min_slow	RQ_plus_slow	RQ_boot_1	RQ_boot_2
    RV  = data(:,1);
    RQ  = data(:,2);
    RJ  = data(:,3);
    BPV = data(:,4);
    RVn = data(:,5);
    RVp = data(:,6);
    
    RVw = mean([RV,mlag(RV,4,mean(RV))],2);
    
    T = length(RV);
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
        
        RV1  = RV1(22:end);
        RV1d = RV1d(22:end);
        RV1w = RV1w(22:end);
        RV1m = RV1m(22:end);
        
        minRV = min(RV1w);  % will use these values in the "insanity filter" below
        maxRV = max(RV1w);
        meanRV = mean(RV1w);
        
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
        
        RV1  = RV1(22:end);
        RV1d = RV1d(22:end);
        RV1w = RV1w(22:end);
        RV1m = RV1m(22:end);
        
        minRV = min(RV1w);  % will use these values in the "insanity filter" below
        maxRV = max(RV1w);
        meanRV = mean(RV1w);
        
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
    
    table7 = nan(4,9);
    for mm=1:9
        table7(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
        table7(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
        
        table7(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
        table7(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window
    end
    table7s(:,:,ss) = table7;
    
end % takes about 13 mins

table7m = mean(table7s,3);

clear info;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F','HARQ-h');
info.rnames = char('Table 7','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
info.width = 120;
if 0
    mprint(table7m,info)
end

table7mm = nan(4,9,2);
for ii=1:9;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table7s(:,ii,ss)./table7s(:,2,ss);
    end
    table7mm(:,ii,1) =   mean(temp,2);
    table7mm(:,ii,2) = median(temp,2);
end
table77mm = nan(8,9);
table77mm(1:2:end,:) = table7mm(:,:,1);
table77mm(2:2:end,:) = table7mm(:,:,2);

info.rnames = char('Table 7','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table77mm,info)

if 0  % printing out more details on the individual stocks
    clear info;
    info.width = 120;
    info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F','HARQ-h');
    info.rnames = char('*MSE-RW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table7s(1,:,:))',info)
    info.rnames = char('*MSE-IW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table7s(2,:,:))',info)
    info.rnames = char('*QLIKE-RW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table7s(3,:,:))',info)
    info.rnames = char('*QLIKE-IW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table7s(4,:,:))',info)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 8: OOS Monthly estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table8s = nan(4,9,num_stocks);
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);  % loading the data for the ss^th stock
    % data is:  RV	BPV	MedRV	TrRV	RV_min	RV_plus	RQ	TPQ	MedRQ	TrRQ	RQ_min	RQ_plus	J	RQ_slow	TPQ_slow	MedRQ_slow	TrRQ_slow	RQ_min_slow	RQ_plus_slow	RQ_boot_1	RQ_boot_2
    RV  = data(:,1);
    RQ  = data(:,2);
    RJ  = data(:,3);
    BPV = data(:,4);
    RVn = data(:,5);
    RVp = data(:,6);
    
    RVm = mean([RV,mlag(RV,21,mean(RV))],2);
    T = length(RV);
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
        
        RV1  = RV1(22:end);
        RV1d = RV1d(22:end);
        RV1w = RV1w(22:end);
        RV1m = RV1m(22:end);
        
        minRV = min(RV1m);  % will use these values in the "insanity filter" below
        maxRV = max(RV1m);
        meanRV = mean(RV1m);
        
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
        
        RV1 = RV1(22:end);
        RV1d = RV1d(22:end);
        RV1w = RV1w(22:end);
        RV1m = RV1m(22:end);
        
        minRV = min(RV1m);  % will use these values in the "insanity filter" below
        maxRV = max(RV1m);
        meanRV = mean(RV1m);
        
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
    
    table8 = nan(4,9);
    for mm=1:9
        table8(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
        table8(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
        
        table8(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
        table8(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window
    end
    table8s(:,:,ss) = table8;
end
table8m = mean(table8s,3);

clear info;
info.width = 120;
info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F','HARQ-h');
info.rnames = char('Table 8','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
if 0
    mprint(table8m,info)
end

table8mm = nan(4,9,2);
for ii=1:9;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table8s(:,ii,ss)./table8s(:,2,ss);
    end
    table8mm(:,ii,1) =   mean(temp,2);
    table8mm(:,ii,2) = median(temp,2);
end
table88mm = nan(8,9);
table88mm(1:2:end,:) = table8mm(:,:,1);
table88mm(2:2:end,:) = table8mm(:,:,2);

info.rnames = char('Table 8','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table88mm,info)

if 0  % printing some more details on the individual stocks
    clear info;
    info.width = 120;
    info.cnames = char('AR','HAR','HAR-J','CHAR','SHAR','ARQ','HARQ','HARQ-F', 'HARQ-h');
    info.rnames = char('*MSE-RW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table8s(1,:,:))',info)
    info.rnames = char('*MSE-IW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table8s(2,:,:))',info)
    info.rnames = char('*QLIKE-RW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table8s(3,:,:))',info)
    info.rnames = char('*QLIKE-IW*','AXP','BA','CAT','CSCO','CVX','DD','DIS','GE','HD','IBM','INTC','JNJ','JPM','KO','MCD','MMM','MRK','MSFT','NKE','PFE','PG','TRV','UNH','UTX','VZ','WMT','XOM');
    mprint(squeeze(table8s(4,:,:))',info)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 9: OOS (daily) comparisons with other RV measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table9s = nan(4,6,num_stocks);
format short g
for ss = 1:num_stocks
    temp_str = ['data1min = xlsread(''',data_dir,'\Data_1min\',dt1min(ss).name,''',''1'',','''b2:g4203'');'];
    evalin('base',temp_str);
    
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
    table9s(:,:,ss) = table9;
    
end
table9m = mean(table9s,3);

clear info;
info.width = 120;
info.cnames = char('HARQ','RV','SS-RV','TS-RV','RK','PA-RV');
info.rnames = char('Table 9','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
if 0
    mprint(table9m,info)
end

table9mm = nan(4,6,2);
for ii=1:6;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table9s(:,ii,ss)./table9s(:,1,ss);
    end
    table9mm(:,ii,1) =   mean(temp,2);
    table9mm(:,ii,2) = median(temp,2);
end
table99mm = nan(8,6);
table99mm(1:2:end,:) = table9mm(:,:,1);
table99mm(2:2:end,:) = table9mm(:,:,2);

info.rnames = char('Table 9','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table99mm,info)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 10: HAR vs HARQ using different RV measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table10s = nan(4,5,num_stocks);
for ss = 1:num_stocks
    temp_str = ['data1min = xlsread(''',data_dir,'\Data_1min\',dt1min(ss).name,''',''1'',','''b2:g4203'');'];
    evalin('base',temp_str);
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
        
        minRV   = min(RV1d);
        maxRV = max(RV1d);
        meanRV = mean(RV1d);    % will use these values in the "insanity filter" below
        
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
        
        minRV   = min(RV1(tt-1000:tt-1));
        maxRV = max(RV1(tt-1000:tt-1));
        meanRV = mean(RV1(tt-1000:tt-1));    % will use these values in the "insanity filter" below
        
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
    for mm=1:10
        table10a(1,mm) = mean( (yhatALL(:,1,1)-yhatALL(:,mm+1,1)).^2 );  % MSE - rolling window
        table10a(2,mm) = mean( (yhatALL(:,1,2)-yhatALL(:,mm+1,2)).^2 );  % MSE - expanding window
        table10a(3,mm) = mean( yhatALL(:,1,1)./yhatALL(:,mm+1,1) - log(yhatALL(:,1,1)./yhatALL(:,mm+1,1)) -1 );  % QLIKE - rolling window
        table10a(4,mm) = mean( yhatALL(:,1,2)./yhatALL(:,mm+1,2) - log(yhatALL(:,1,2)./yhatALL(:,mm+1,2)) -1 );  % QLIKE - expanding window
    end
    table10=table10a(:,2:2:end)./table10a(:,1:2:end);
    table10s(:,:,ss) = table10;
end

table10m = mean(table10s,3);

clear info;
info.cnames = char('RV','SS-RV','TS-RV','RK','PA-RV');
info.rnames = char('Table 10','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
if 0 
    mprint(table10m,info)
end


table10mm = nan(4,5,2);

for ii=1:5;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table10s(:,ii,ss);
    end
    table10mm(:,ii,1) =   mean(temp,2);
    table10mm(:,ii,2) = median(temp,2);
end
table1010mm = nan(8,5);
table1010mm(1:2:end,:) = table10mm(:,:,1);
table1010mm(2:2:end,:) = table10mm(:,:,2);

info.rnames = char('Table 10','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table1010mm,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 11: Different IQ estimators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table11s = nan(4,6,num_stocks);

% Tables 11 and 12 use the SP500 index data again.
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);  % loading the data for the ss^th stock
    
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
    table11s(:,:,ss)=table11;
end
table11m = mean(table11s,3);

clear info;
info.cnames = char('RQ','TPQ','MedRQ','TrRQ','RQ15min','Boot');
info.rnames = char('Table 11','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
if 0
    mprint(table11m,info)
end

table11mm = nan(4,6,2);
for ii=1:6;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table11s(:,ii,ss)./table11s(:,1,ss);
    end
    table11mm(:,ii,1) =   mean(temp,2);
    table11mm(:,ii,2) = median(temp,2);
end
table1111mm = nan(8,6);
table1111mm(1:2:end,:) = table11mm(:,:,1);
table1111mm(2:2:end,:) = table11mm(:,:,2);

info.rnames = char('Table 11','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table1111mm,info)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 12: Different HARQ specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table12s = nan(4,7,num_stocks);
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);  % loading the data for the ss^th stock
    
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
    
    T = size(RV,1);
    
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
    table12s(:,:,ss) = table12;
end
table12m = mean(table12s,3);

clear info;
info.width = 120;
info.cnames = char('RQ','RQ^1/2','RQ^-1/2','RQ^-1','Log RQ','HAR','HARQ');
info.rnames = char('Table 12','MSE-RW','MSE-IW','QLIKE-RW','QLIKE-IW');
info.fmt    = '%10.4f';
sprintf('\n\n');
if 0
    mprint(table12m,info)
end

table12mm = nan(4,7,2);
for ii=1:7;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table12s(:,ii,ss)./table12s(:,2,ss);
    end
    table12mm(:,ii,1) =   mean(temp,2);
    table12mm(:,ii,2) = median(temp,2);
end
table1212mm = nan(8,7);
table1212mm(1:2:end,:) = table12mm(:,:,1);
table1212mm(2:2:end,:) = table12mm(:,:,2);

info.rnames = char('Table 12','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table1212mm,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 13: Alternative Q-model in sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table13 = zeros(3,6);

for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);
    
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
    T = length(RV);
    
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
    
    
    hhatALLis = [hhatHARJ,hhatHARQJ,hhatCHAR,hhatCHARQ,hhatSHAR,hhatSHARQ];
    for ii=1:size(hhatALLis,2);
        table13(1,ii) = table13(1,ii) + (1-cov((RV(23:end)-hhatALLis(:,ii)).*(hhatALLis(:,ii)>=0))/cov(RV(23:end)))/num_stocks;            % r2
        table13(2,ii) = table13(2,ii) + ( mean( ((RV(23:end)-hhatALLis(:,ii)).^2).*(hhatALLis(:,ii)>=0)))/num_stocks;                              % MSE
        table13(3,ii) = table13(3,ii) + ((mean( (RV(23:end)./hhatALLis(:,ii) - log(RV(23:end)./hhatALLis(:,ii)) -1 ).*(hhatALLis(:,ii)>=0) ) ))/num_stocks;   % QLIKE
    end
end

clear info;
info.cnames = char('HAR-J', 'HARQ-J', 'CHAR', 'CHARQ', 'SHAR', 'SHARQ');
info.rnames = char('Table 13','R2','MSE','QLIKE');
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
table14s = nan(4,4,num_stocks);
for ss = 1:num_stocks
    temp_str = ['data = xlsread(''',data_dir,'\Data_5min\',dt(ss).name,''',''1'',','''b2:l4203'');'];
    evalin('base',temp_str);  % loading the data for the ss^th stock
    
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
    
    T = size(RV,1);
    
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
        tempSHARQ = ols(RV1(2:end),[ones(window-22,1),RV1w(1:end-1),RV1m(1:end-1),RVp1(1:end-1),RVn1(1:end-1),RVp1(1:end-1).*sqrtRQd(1:end-1), RVn1(1:end-1).*sqrtRQd(1:end-1)]);
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
    table14s(:,:,ss) = table14;
end

table14mm = nan(4,4,2);
for ii=1:4;
    temp = nan(4,num_stocks);
    for ss=1:num_stocks
        temp(:,ss) = table14s(:,ii,ss);
    end
    table14mm(:,ii,1) = mean(temp,2);
    table14mm(:,ii,2) = median(temp,2);
end
table1414mm = nan(8,4);
table1414mm(1:2:end,:) = table14mm(:,:,1);
table1414mm(2:2:end,:) = table14mm(:,:,2);

info.cnames = char('HARQ', 'HARQ-J', 'CHARQ', 'SHARQ');
info.rnames = char('Table 14','MSE-RW-Avg','MSE-RW-Med','MSE-IW-Avg','MSE-IW-Med','QLIKE-RW-Avg','QLIKE-RW-Med','QLIKE-IW-Avg','QLIKE-IW-Med');
mprint(table1414mm,info)


toc/60