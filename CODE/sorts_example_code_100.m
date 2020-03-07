% Matlab code to re-produce some of the results in:
%
% "Monotonicity in Asset Returns: New Tests with Applications to the Term
% Structure, the CAPM and Portfolio Sorts"
%
% Paper available at:
% http://www.economics.ox.ac.uk/members/andrew.patton/research.html
%
%
%  Andrew Patton and Allan Timmermann
%
%  28 April 2009


% This paper is a revised version of the paper previously circulated as:
% 
%  "Portfolio Sorts and Tests of Cross-Sectional Patterns in Expected Returns"
%    by Andrew J. Patton and Allan Timmermann


bootreps = 1000
% takes about 6.7 minutes to run on a 2.4GHz machine with 3GB of memory, with bootreps=1000, and wolak_reps=0

wolak_sims = 100
% this should be set to at least 100, and preferably 1000, if the Wolak tests are to be used. (if they are to be ignored then set to 2, say.)
% takes  18.0 minutes for wolak_sims=100
% takes 169.4 minutes for wolak_sims=1000
% in our application the results did not change much at all for sims=100 or sims=1000.


tic;

data_directory = 'c:\core\work\sorts\';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLICATION 1: CAPM BETA SORTED PORTFOLIOS - average returns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% all stocks on NYSE, NASDAQ and AMEX
evalin('base',['load ''',data_directory,'bt_a_10_vw.txt'' -ascii;']);
data1 = bt_a_10_vw(:,2:end)*100;
dates = bt_a_10_vw(:,1);
[dates(1),dates(end)]  %  196307      200112
[T,n] = size(data1)  % 462  10

% plotting the data
figure(1),plot((1:10)',mean(data1),'o-','LineWidth',1);hold on
for ii=2:10;
    if mean(data1(:,ii))<mean(data1(:,ii-1));
        plot([ii-1,ii],mean(data1(:,ii-1:ii)),'r.-');
    end
end
set(gca,'XTickLabel',{'Low';'2';'3';'4';'5';'6';'7';'8';'9';'High'});
axis([0.8,10.2,0.39,0.61]),...
    ylabel('Average return'),...
    title('Value-weighted past beta portfolio returns, 1963-2001');hold off;


% running the tests of monotonicity
direction=1;
block_length=10;
rand_state=1234;  % note that if a different rand_state value is chosen then slightly different values will be returned.
table1b = -999.99*ones(1,9);
data123 = data1;

% top minus bottom
table1b(1,1) = mean(data123(:,end)-data123(:,1));

% t-tests and MR tests
temp = MR_test_22(data123,bootreps,direction,block_length,rand_state);
table1b(1,2:5) = temp(:,2)';

% "up" and "down" tests
temp = MR_up_down_test_20(data123,bootreps,block_length,rand_state);
table1b(1,6:7)= temp(3:4,2);  % using the U/D test on abs values, studentised

% wolak test
temp = wolak_test_3((data123(:,2:end)-data123(:,1:end-1))*direction,wolak_sims);
table1b(1,8) = temp(1);

% Bonferroni test
temp = Bonferroni_MR_test1((data123(:,2:end)-data123(:,1:end-1))*direction);
table1b(1,9) = min(temp(1),1);  % the way I report "bonferroni p-values" allows these to be greater than 1, so set max at 1

% creating table 2A
table1a = -999.99*ones(2,10);
table1a(1,:) = mean(data1);
table1a(2,:) = std(data1);
clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('Low','2','3','4','5','6','7','8','9','High');
info.rnames = strvcat('.','Mean','Std Dev');
sprintf(['Table 2A: Average returns on CAPM beta decile portfolios'])
mprint(table1a,info)

% printing the results of the monotonicity tests
clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('top-bottom','t-stat','t-pval','MR-pval','MRall-pval','UP-pval','DOWN-pval','Wolak-pval','Bonf-pval');
sprintf(['Table 2B: Test of monotonicity for returns on CAPM beta portfolios'])
mprint(table1b,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLICATION 2: CAPM BETA SORTED PORTFOLIOS - ex-post betas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loading in the market and rf rate
data2 = xlsread([data_directory,'F-F_Research_Data_Factors.xls'],'Clean','A446:E907');
mkt = data2(:,2);
rf = data2(:,5);

% getting the ex-post betas on the (ex-ante) CAPM beta decile portfolios
table1c = -999.99*ones(1,10);
for ii=1:10;
    temp = ols(data1(:,ii),[ones(T,1),mkt]);
    table1c(ii) = temp.beta(2);
end

% printing as a table
clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('Low','2','3','4','5','6','7','8','9','High');
sprintf(['Table 2C: Estimated ex-post betas for (ex-ante) CAPM beta portfolios'])
mprint(table1c,info)

% running the tests of monotonicity on the estimated betas
direction=1;
block_length=10;
rand_state=1234;  % note that if a different rand_state value is chosen then slightly different values will be returned.
temp = MR_test_SLOPE_2(data1,[mkt,ones(T,1)],0,bootreps,direction,block_length,rand_state);
table1d = temp(:,2)';

clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('top-bottom','t-stat','t-pval','MR-pval','MRall-pval');
sprintf(['Table 2D: Tests of monotonicity on ex-post betas of (ex-ante) CAPM beta portfolios'])
mprint(table1d,info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLICATION 1: TERM PREMIUM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1 = xlsread([data_directory,'fama_files.xls'],'filled','A2:M463');
dates = data1(:,1);
data1 = data1(:,2:end);
[dates(1),dates(end)]  % 19630731    20011231
[T,k] = size(data1)  %  462    12, goes from 1 month to 12 month maturities (first col is short maturity, last col is long maturity)

% pulling out the sub-sample data
[dates(7),dates(114)]  %   19640131    19721229, these are the dates of the first and last obs of the first sub-sample in BRSW
termprem0 = 100*(data1(7:end,2:end-1)-data1(7:end,1)*ones(1,k-2));  % BRSW drop the 12-month maturity bond as there are many missing values (same in our data set)
termprem1 = 100*(data1(7:114,2:end-1)-data1(7:114,1)*ones(1,k-2));
termprem2 = 100*(data1(115:end,2:end-1)-data1(115:end,1)*ones(1,k-2));

% plotting the average term premia over the full sample
figure(2),plot((2:11)',mean(termprem0),'o-','LineWidth',1);hold on
for ii=2:10;
    if mean(termprem0(:,ii))<mean(termprem0(:,ii-1));
        plot([ii,ii+1],mean(termprem0(:,ii-1:ii))','r.-');
    end
end
title('US T-bill term premia, 1964-2001'),...
    axis([1.5,11.5,0,0.1]),...
    xlabel('Maturity (months)'),ylabel('Average term premium');hold off;

% running the tests of monotonicity
direction=1;
block_length=10;
rand_state=1234;  % note that if a different rand_state value is chosen then slightly different values will be returned. 
table2b = -999.99*ones(3,9);
for ii=1:3;
    if ii==1;
        data123 = termprem0;
    elseif ii==2
        data123 = termprem1;
    else
        data123 = termprem2;
    end

    % top minus bottom
    table2b(ii,1) = mean(data123(:,end)-data123(:,1));

    % t-tests and MR tests
    temp = MR_test_22(data123,bootreps,direction,block_length,rand_state);
    table2b(ii,2:5) = temp(:,2)';
    
    % "up" and "down" tests
    temp = MR_up_down_test_20(data123,bootreps,block_length,rand_state);
    table2b(ii,6:7)= temp(3:4,2);  % using the U/D test on abs values, studentised

    % wolak test
    temp = wolak_test_3((data123(:,2:end)-data123(:,1:end-1))*direction,wolak_sims);
    table2b(ii,8) = temp(1);
    
    % Bonferroni test
    temp = Bonferroni_MR_test1((data123(:,2:end)-data123(:,1:end-1))*direction);
    table2b(ii,9) = min(temp(1),1);  % the way I report "bonferroni p-values" allows these to be greater than 1, so set max at 1
end

% computing average term premia in the three samples
table2a = -999.99*ones(3,10);
table2a(1,:) = mean(termprem0);
table2a(2,:) = mean(termprem1);
table2a(3,:) = mean(termprem2);

% this is Table 1 in the paper
% printing table in nice format
clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('2','3','4','5','6','7','8','9','10','11');
info.rnames = strvcat('Sample period','1964-2001','1964-1972','1973-2001')
sprintf(['Table 2A: Average term premia'])
mprint(table2a,info)
info.cnames = strvcat('top-bottom','t-stat','t-pval','MR-pval','MRall-pval','UP-pval','DOWN-pval','Wolak-pval','Bonf-pval');
sprintf(['Table 2B: Test statistics for monthly term premia (Jan 1964- December 2001)'])
mprint(table2b,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLICATION 4: ONE-WAY PORTFOLIO SORTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data originally downloaded from Ken French's web site

% MONTHLY RETURNS, 1963.07 to 2006.12 (522 months)
data1 = xlsread([data_directory,'PT08_data_ALL.xls'],'ME',      'K446:T967');
data2 = xlsread([data_directory,'PT08_data_ALL.xls'],'BE-ME',   'K446:T967');
data3 = xlsread([data_directory,'PT08_data_ALL.xls'],'CF-P',    'K146:T667');
data4 = xlsread([data_directory,'PT08_data_ALL.xls'],'E-P',     'K146:T667');
data5 = xlsread([data_directory,'PT08_data_ALL.xls'],'D-P',     'K434:T955');
data6 = xlsread([data_directory,'PT08_data_ALL.xls'],'Momentum','B440:K961');
data7 = xlsread([data_directory,'PT08_data_ALL.xls'],'STR',     'B451:K972');
data8 = xlsread([data_directory,'PT08_data_ALL.xls'],'LTR',     'B392:K913');
dataALL = [data1,data2,data3,data4,data5,data6,data7,data8];

% MONTHLY RETURNS, all available data (will be different lengths)
data1b = xlsread([data_directory,'PT08_data_ALL.xls'],'ME',      'K2:T967');
data2b = xlsread([data_directory,'PT08_data_ALL.xls'],'BE-ME',   'K2:T967');
data3b = xlsread([data_directory,'PT08_data_ALL.xls'],'CF-P',    'K2:T667');
data4b = xlsread([data_directory,'PT08_data_ALL.xls'],'E-P',     'K2:T667');
data5b = xlsread([data_directory,'PT08_data_ALL.xls'],'D-P',     'K2:T955');
data6b = xlsread([data_directory,'PT08_data_ALL.xls'],'Momentum','B2:K961');
data7b = xlsread([data_directory,'PT08_data_ALL.xls'],'STR',     'B2:K972');
data8b = xlsread([data_directory,'PT08_data_ALL.xls'],'LTR',     'B2:K913');


% creating Table 3a, using the post 1963 data
table3a = -999.99*ones(11,8);
for ii=1:8;
     table3a(1:10,ii) = mean(dataALL(:,(ii-1)*10+1:ii*10))';
     table3a(11,ii) = table3a(10,ii)-table3a(1,ii);
end

clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('ME','BE-ME','CF-P','E-P','D-P','Momentum','ST reversal','LT reversal');
info.rnames = strvcat('.','Low','2','3','4','5','6','7','8','9','High','High-Low');
sprintf(['Table 3A: Average returns, 1963-2006']);
mprint(table3a,info)

% running the tests of monotonicity on the post-1963 portfolio data
block_length=10;
rand_state=1234;  % note that if a different rand_state value is chosen then slightly different values will be returned.
direction=[-1;1;1;1;1;1;-1;-1];  % directions are different depending on the portfolio sort we're looking at


table3b = -999.99*ones(9,8);
for ii=1:8;
    data123 = dataALL(:,(ii-1)*10+1:ii*10);

    % top minus bottom
    table3b(1,ii) = mean(data123(:,end)-data123(:,1));

    % t-tests and MR tests
    temp = MR_test_22(data123,bootreps,direction(ii),block_length,rand_state);
    table3b(2:5,ii) = temp(:,2)';
    
    % "up" and "down" tests
    temp = MR_up_down_test_20(data123,bootreps,block_length,rand_state);
    table3b(6:7,ii)= temp(3:4,2);  % using the U/D test on abs values, studentised

    % wolak test
    temp = wolak_test_3((data123(:,2:end)-data123(:,1:end-1))*direction(ii),wolak_sims);
    table3b(8,ii) = temp(1);
    
    % Bonferroni test
    temp = Bonferroni_MR_test1((data123(:,2:end)-data123(:,1:end-1))*direction(ii));
    table3b(9,ii) = min(temp(1),1);  % the way I report "bonferroni p-values" allows these to be greater than 1, so set max at 1
end

% printing table 3b
clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('ME','BE-ME','CF-P','E-P','D-P','Momentum','ST reversal','LT reversal');
info.rnames = strvcat('.','top-bottom','t-stat','t-pval','MR-pval','MRall-pval','UP-pval','DOWN-pval','Wolak-pval','Bonf-pval');
sprintf(['Table 3B: Tests of monotonicity, 1963-2006'])
mprint(table3b,info)


% creating Table 3C, using the FULL SAMPLE data
table3a = -999.99*ones(11,8);
for ii=1:8;
    if     ii==1; data123 = data1b;
    elseif ii==2; data123 = data2b;
    elseif ii==3; data123 = data3b;
    elseif ii==4; data123 = data4b;
    elseif ii==5; data123 = data5b;
    elseif ii==6; data123 = data6b;
    elseif ii==7; data123 = data7b;
    elseif ii==8; data123 = data8b; end
    table3c(1:10,ii) = mean(data123)';
    table3c(11,ii) = table3c(10,ii)-table3c(1,ii);
end

clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('ME','BE-ME','CF-P','E-P','D-P','Momentum','ST reversal','LT reversal');
info.rnames = strvcat('.','Low','2','3','4','5','6','7','8','9','High','High-Low');
sprintf(['Table 3C: Average returns, full sample'])
mprint(table3c,info)

% running the tests of monotonicity
table3d = -999.99*ones(9,8);
for ii=1:8;
    if     ii==1; data123 = data1b;
    elseif ii==2; data123 = data2b;
    elseif ii==3; data123 = data3b;
    elseif ii==4; data123 = data4b;
    elseif ii==5; data123 = data5b;
    elseif ii==6; data123 = data6b;
    elseif ii==7; data123 = data7b;
    elseif ii==8; data123 = data8b; end

    % top minus bottom
    table3d(1,ii) = mean(data123(:,end)-data123(:,1));

    % t-tests and MR tests
    temp = MR_test_22(data123,bootreps,direction(ii),block_length,rand_state);
    table3d(2:5,ii) = temp(:,2)';
    
    % "up" and "down" tests
    temp = MR_up_down_test_20(data123,bootreps,block_length,rand_state);
    table3d(6:7,ii)= temp(3:4,2);  % using the U/D test on abs values, studentised

    % wolak test
    temp = wolak_test_3((data123(:,2:end)-data123(:,1:end-1))*direction(ii),wolak_sims);
    table3d(8,ii) = temp(1);
    
    % Bonferroni test
    temp = Bonferroni_MR_test1((data123(:,2:end)-data123(:,1:end-1))*direction(ii));
    table3d(9,ii) = min(temp(1),1);  % the way I report "bonferroni p-values" allows these to be greater than 1, so set max at 1
end

% printing table 3d
clear info;
info.fmt = '%10.3f';
info.cnames = strvcat('ME','BE-ME','CF-P','E-P','D-P','Momentum','ST reversal','LT reversal');
info.rnames = strvcat('.','top-bottom','t-stat','t-pval','MR-pval','MRall-pval','UP-pval','DOWN-pval','Wolak-pval','Bonf-pval');
sprintf(['Table 3D: Tests of monotonicity, full sample'])
mprint(table3d,info)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLICATION 5: TWO-WAY PORTFOLIO SORTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table4a = -999.99*ones(7,7);
table4b = -999.99*ones(7,7);

% MONTHLY RETURNS, using post-1963 data
% data originally downloaded from Ken French's web site
data12b = xlsread([data_directory,'PT08_data_ALL.xls'],'ME_BEME',  'B446:Z967');
data16b = xlsread([data_directory,'PT08_data_ALL.xls'],'ME_MOM',   'B440:Z961');

% there are some missing values for some of these (early in the sample) so will drop those months
data12b = data12b(find(min(data12b,[],2)>-99.99),:);
data16b = data16b(find(min(data16b,[],2)>-99.99),:);

% creating a 3D matrix of returns on the two-way sorted portfolios, which makes it easier to see what's going on
data12 = nines(5,5,size(data12b,1));
data16 = nines(5,5,size(data16b,1));
for ii=1:5;
    data12(ii,:,:) = data12b(:,5*(ii-1)+1:5*ii)';
    data16(ii,:,:) = data16b(:,5*(ii-1)+1:5*ii)';
end

% computing the mean returns on these portfolios
table4a(1:5,1:5) = mean(data12,3);
table4b(1:5,1:5) = mean(data16,3);

% running the monotonicity tests

data12r = data12(:,end:-1:1,:);  % reversing order of columns 
data16r = data16(:,end:-1:1,:);  % reversing order of columns
% this is so we are looking for a pattern with expected returns decreasing from top to bottom, and left to right

% need to create a matrix with all the differences to be tested
data12a = []; data16a = []; 
for ii=1:5;
    data12a = [data12a,squeeze(data12r(ii,2:end,:)-data12r(ii,1:end-1,:))'];  % all the column differences
    data12a = [data12a,squeeze(data12r(2:end,ii,:)-data12r(1:end-1,ii,:))'];  % all the row differences

    data16a = [data16a,squeeze(data16r(ii,2:end,:)-data16r(ii,1:end-1,:))'];  % all the column differences
    data16a = [data16a,squeeze(data16r(2:end,ii,:)-data16r(1:end-1,ii,:))'];  % all the row differences
end

% overall test for monotonicity
temp1 = MR_test_22(data12a,bootreps,-2,block_length,rand_state);
temp2 = MR_test_22(data16a,bootreps,-2,block_length,rand_state);

table4a(7,7) = temp1(3,2);  % using studentised MR test
table4b(7,7) = temp2(3,2);

% now doing the tests for individual rows and columns
for ii=1:5;
    temp = MR_test_22(squeeze(data12r(2:end,ii,:)-data12r(1:end-1,ii,:))',bootreps,-2,block_length,rand_state);
    table4a(6,6-ii) = temp(3,2);  % i re-order the column that the p-value is placed into to undo the re-ordering i used above
    temp = MR_test_22(squeeze(data12r(ii,2:end,:)-data12r(ii,1:end-1,:))',bootreps,-2,block_length,rand_state);
    table4a(ii,6) = temp(3,2);
    
    temp = MR_test_22(squeeze(data16r(2:end,ii,:)-data16r(1:end-1,ii,:))',bootreps,-2,block_length,rand_state);
    table4b(6,6-ii) = temp(3,2);
    temp = MR_test_22(squeeze(data16r(ii,2:end,:)-data16r(ii,1:end-1,:))',bootreps,-2,block_length,rand_state);
    table4b(ii,6) = temp(3,2);
end

% now doing joint tests for all rows, or all columns
temp2 = [];
temp3 = [];
for ii=1:5;
    temp2 = [temp2,squeeze(data12r(2:end,ii,:)-data12r(1:end-1,ii,:))'];
    temp3 = [temp3,squeeze(data12r(ii,2:end,:)-data12r(ii,1:end-1,:))'];
end
temp = MR_test_22(temp2,bootreps,-2,block_length,rand_state);
table4a(7,3) = temp(3,2);
temp = MR_test_22(temp3,bootreps,-2,block_length,rand_state);
table4a(3,7) = temp(3,2);
temp2 = [];
temp3 = [];
for ii=1:5;
    temp2 = [temp2,squeeze(data16r(2:end,ii,:)-data16r(1:end-1,ii,:))'];
    temp3 = [temp3,squeeze(data16r(ii,2:end,:)-data16r(ii,1:end-1,:))'];
end
temp = MR_test_22(temp2,bootreps,-2,block_length,rand_state);
table4b(7,3) = temp(3,2);
temp = MR_test_22(temp3,bootreps,-2,block_length,rand_state);
table4b(3,7) = temp(3,2);


clear info;
info.fmt = '%10.3f';
info.rnames = strvcat('.','Small','2','3','4','Big','MR pval','Joint MR pval');
info.cnames = strvcat('Growth','2','3','4','Value','MR pval','Joint MR pval');
sprintf(['Table 4A: Conditional and joint monotonicity tests for double-sorted portfolios, ME x BE/ME'])
mprint(table4a,info)

info.cnames = strvcat('Losers','2','3','4','Winners','MR pval','Joint MR pval');
sprintf(['Table 4B: Conditional and joint monotonicity tests for double-sorted portfolios, ME x M''tum'])
mprint(table4b,info)

toc