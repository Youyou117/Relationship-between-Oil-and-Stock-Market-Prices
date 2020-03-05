%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADING IN THE DATA
[num_data txt_data]= xlsread('30MIN数据.xlsx');   

date = txt_data(2:end,1);


datet=datetime(date);
rets1(:,1:2)=num_data(:,1:2);
T = length(rets1)  % 2589

sectors = txt_data(1,2:end);  %%%各金融机构的名称

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC SUMMARY STATISTICS AND FIGURES
%%%说明，在进行jbtest、adftest和archtest时，需移除路径D:\MATLAB\R2016b\toolbox\jplv7\Ucsd_garch\Garch
%%%也就是移除jplv7\Ucsd_garch\Garch

table1a = nan(size(rets1,2),11);

table1aP=nan(size(rets1,2),11);   %%对应P值
table1a(:,1) = mean(rets1)';
table1a(:,2) = max(rets1)';
table1a(:,3) = min(rets1)';
table1a(:,4) = std(rets1)';
table1a(:,5) = skewness(rets1)';
table1a(:,6) = kurtosis(rets1)';
for i=1:size(rets1,2)
    [h jbp jbstat]=jbtest(rets1(:,i));
    table1a(i,7) =jbstat;
    table1aP(i,7) =jbp;
    [h adfp adfstat]=adftest(rets1(:,i));
    table1a(i,8) =adfstat;
    table1aP(i,8) =adfp;
    [h,pValue,lbstat,cValue] = lbqtest(rets1(:,i),'lags',20);  %% Ljung-Box statistics of returns,Q(20)
    table1a(i,9) =lbstat;
    table1aP(i,9) =pValue;
    [h,pValue,lb2stat,cValue] = lbqtest(rets1(:,i).^2,'lags',20);  %% Ljung-Box statistics of returns,Q^2(20)
    table1a(i,10) =lb2stat;
    table1aP(i,10) =pValue;
    [h,pValue,arstat,cValue] = archtest(rets1(:,i),'Lags',20);  %% archtest,Engle's LM test for heteroskedasticity computed using 20 lags
    table1a(i,11) =arstat;
    table1aP(i,11) =pValue;
end



info.fmt = '%10.3f';
txt_data(1,2:end)
info.rnames = strvcat('.','石油期货','中国石油');
info.cnames = strvcat('Mean','Max','Min','Std dev','Skewness','Kurtosis','Jarque-Bera','ADF','Q(20)','Q^2(20)','ARCH(20)');
sprintf(['Table 1a: Summary statistics'])
mprint(table1a,info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBTAINING MODELS FOR THE CONDITIONAL MEAN AND VARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First: conditional mean

resids = nan(T,size(rets1,2));
mean_order = nan(size(rets1,2),3);
tic;
for ii=1:size(rets1,2)
    [theta1,sig21,vcv1,order1,resids1] = ARMAX_opt(rets1(:,ii),3,3,'AIC');  % takes about 15 seconds per variable
    mean_order(ii,1:2) = order1';
    mean_order(ii,3) = 1-sig21/cov(rets1(:,ii));
    resids(:,ii) = [zeros(max(order1),1);resids1];
    [ii,toc]
end
toc  % takes 15 seconds for each series
mean_order 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second: the conditional variance
% considerthese models:
% Const_vol, ARCH(1), GARCH(1,1), GJR-GARCH(1,1,1), GARCH(2,2), GJR-GARCH(2,2,2)
 info.fmt = '%10.3f';     

vol_LL = nan(size(rets1,2),7,2); % will use Normal likelihood to compare these vol models
hhat_ALL = nan(T,size(rets1,2),7);
tic;
GARCHparamsn = nan(7,size(rets1,2),7);  %%第一个7表示最多7个参数，用于接收7种模型的参数


for ii=1:size(rets1,2)
    % model 1: constant volatility
    hhat_ALL(:,ii,1) = cov(resids(:,ii))*ones(T,1);
    
    % model 2: ARCH(1)
    [parameters2, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),1,0,0,'GJRGARCH','NORMAL',[],[cov(resids(:,ii))*(1-0.3);0.3]);parameters2
    hhat_ALL(:,ii,2) = hhat1;
    GARCHparamsn(1:2,ii,2)=parameters2;
    % model 3: GARCH(1,1)
    [parameters3, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),1,0,1,'GJRGARCH','NORMAL',[],[cov(resids(:,ii))*0.05;0.05;0.9]);parameters3
    hhat_ALL(:,ii,3) = hhat1;
    GARCHparamsn(1:3,ii,3)=parameters3;
    % model 4: GJR-GARCH(1,1,1)
    [parameters4, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),1,1,1,'GJRGARCH','NORMAL',[],[parameters3(1:2);0;parameters3(3)]);parameters4
    hhat_ALL(:,ii,4) = hhat1;
    
    % saving these parametesr for the KS and CvM tests below, as this model turns out to win according to the BIC
   GARCHparamsn(1:4,ii,4)=parameters4;
    
    % model 5: ARCH(2)
    [parameters5, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),2,0,0,'GJRGARCH','NORMAL',[],[cov(resids(:,ii))*(1-0.7);0.5;0.2]);parameters5
    hhat_ALL(:,ii,5) = hhat1;
    GARCHparamsn(1:3,ii,5)=parameters5;
    % model 6: GARCH(2,2)
    [parameters6, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),2,0,2,'GJRGARCH','NORMAL',[],[parameters3(1:2);0;parameters3(3);0]);parameters6
    hhat_ALL(:,ii,6) = hhat1;
    GARCHparamsn(1:5,ii,6)=parameters6;
    % model 7: GJR-GARCH(2,2,2)
    [parameters7, likelihood1, ~, ~, hhat1] = multigarch_AJP(resids(:,ii),2,2,2,'GJRGARCH','NORMAL',[],[parameters4(1:2);0;parameters4(3);0;parameters4(4);0]);parameters7
    hhat_ALL(:,ii,7) = hhat1;
    GARCHparamsn(1:7,ii,7)=parameters7;
    % now computing the mean log-like of each of these models (doing it here, rather than using output from multigarch so that I am sure that the constant vol model is being compared in the right way)
    vol_params = [1,2,3,4,3,5,7];  % number of params in the above models
    for mm=1:7;
        vol_LL(ii,mm,1) = -1/2*log(2*pi) -1/2*mean(log(hhat_ALL(:,ii,mm))) - 1/2*mean( (resids(:,ii).^2)./hhat_ALL(:,ii,mm) );
        vol_LL(ii,mm,2) = -2*vol_LL(ii,mm,1) + log(T)/T*vol_params(mm);  % BIC
        
    end
    
    [ii,toc]
end

toc  
vol_LL


vol_order = nan(size(rets1,2),3);
hhat_opt = nan(T,size(rets1,2));
orders_all = [[0,0,0];[1,0,0];[1,1,0];[1,1,1];[2,0,0];[2,2,0];[2,2,2]];
for ii=1:size(rets1,2)
    temp = find(vol_LL(ii,:,2)==min(vol_LL(ii,:,2)));
    vol_order(ii,:) = orders_all(temp,:);
    hhat_opt(:,ii) = hhat_ALL(:,ii,temp);  % conditional variance from the BIC-optimal model
end
vol_order
% both choose garch(1,1).
stdresids = resids./sqrt(hhat_opt);

% scatter plots of data and std resids
figure(151),plot([-0.15,0.15],[0,0],'k--',[0,0],[-0.15,0.15],'k--','LineWidth',2);hold on;
plot(rets1(:,1),rets1(:,2),'bo'),grid on;
xlabel('石油期货 return'),ylabel('中国石油 return'),axis([-13.5,13.5,-13.5,13.5]),title('Daily returns on 石油期货 and 中国石油');


figure(152),plot([-15,15],[0,0],'k--',[0,0],[-15,15],'k--','LineWidth',2);hold on;
plot(stdresids(:,1),stdresids(:,2),'bo'),grid on;
xlabel('PAB'),ylabel('BNB'),axis([-8,8,-8,8]),title('Standardized residuals for PAB and BNB');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SKEW T MODELS FOR THE STANDARDIZED RESIDUALS

options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'DiffMaxChange',Inf,'DiffMinChange',0,'Algorithm','active-set');

outSKEWT = nan(size(rets1,2),2);  % params
lower = [2.1, -0.99];
upper = [Inf, 0.99 ];
theta0 = [6;0];
for ii=1:size(rets1,2);
    theta1 = fmincon('skewtdis_LL',theta0,[],[],[],[],lower,upper,[],options,stdresids(:,ii));
    outSKEWT(ii,:) = theta1';
end
outSKEWT

outSKEWTandstd = nan(size(rets1,2)*2,2);  % params and standard deviation
for ii=1:size(rets1,2);
    theta1 = fmincon('skewtdis_LL',theta0,[],[],[],[],lower,upper,[],options,stdresids(:,ii));
    scoresc = LLgrad_1('skewtdis_LLa',theta1,stdresids(:,ii));
    Hc1 = hessian('skewtdis_LL',theta1,stdresids(:,ii))/T;  % used in naive VCV matrix for  SKEW T params
    
    Bc = newey_west(scoresc,floor(T^(1/3)));
    Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
    [theta1, sqrt(diag(Vc2))]';
    outSKEWTandstd((2*ii-1):(2*ii),:) =[theta1, sqrt(diag(Vc2))]';
end


tFreedom= nan(size(rets1,2),1);  % 数据对应t分布的自由度
for ii=1:size(rets1,2);
    [nhat,nci] = mle(stdresids(:,ii),'pdf',@tpdf,'start',10);
    tFreedom(ii)=nhat;
end

tFreedomandstd= nan(size(rets1,2)*2,1);  % 数据对应t分布的自由度及其标准差
for ii=1:size(rets1,2);
    [nhat,nci] = mle(stdresids(:,ii),'pdf',@tpdf,'start',10);
    tFreedomandstd(2*ii-1)=nhat;
    tFreedomandstd(2*ii)=(nci(2)-nhat)/1.96;
end

Uedf = empiricalCDF(stdresids);  % prob integral transforms using the empirical cdf
Uskewt = nan(T,size(rets1,2));
for ii=1:size(rets1,2);
    Uskewt(:,ii) = skewtdis_cdf(stdresids(:,ii),outSKEWT(ii,1),outSKEWT(ii,2));
    Ustudent(:,ii)=tcdf(stdresids(:,ii),tFreedom(ii));
end
Uall = Uedf;
Uall(:,:,2) = Uskewt;    % so Uall contains nonparam U's first, then SKEWt U's
Uall(:,:,3) = Ustudent;  %%%tudent's T cumulative distribution function 

size(Uall)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%arma-garch模型估计


%%%%%arma模型系数估计
%%%%%说明，需移除路径D:\MATLAB\R2016b\toolbox\jplv7\Ucsd_garch\Garch
info.fmt = '%10.3f';
pqcoef=nan(size(rets1,2)*2,7+3+2+1+4);  %%第一列为截距项，第2:4为ar系数，第5:7为ma系数,8:10为garch估计，第11列为偏t分布的自由度，12列为对称性参数
                                        %%%第13列为logL，14:17列为其它检验
tic
for i=1:size(rets1,2)
    pq=mean_order(i,1:2);   %%分别为ar的阶数，ma的阶数,garch（1,1）
    if pq(1)==0 & pq(2)==0  %%ar
        Mdl = arima('Variance',garch(1,1));
        [EstMdl,EstParamCov,logL,info]  = estimate(Mdl,rets1(:,i));
        pqcoef(2*i-1,1)=info.X(1);
        pqcoef(2*i,1)=sqrt(EstParamCov(1,1));
        pqcoef(2*i-1,8:10)=info.X(2:4);
        pqcoef(2*i,8)=sqrt(EstParamCov(2,2));
        pqcoef(2*i,9)=sqrt(EstParamCov(3,3));
        pqcoef(2*i,10)=sqrt(EstParamCov(4,4));
        pqcoef((2*i-1):(2*i),11)= outSKEWTandstd((2*i-1):(2*i),1);
        pqcoef((2*i-1):(2*i),12)= outSKEWTandstd((2*i-1):(2*i),2);
        
        pqcoef(2*i-1,13)=-logL;
         
    elseif pq(1)>0 & pq(2)==0
        Mdl = arima('ARLags',[1:pq(1)],'Variance',garch(1,1));
        [EstMdl,EstParamCov,logL,info]  = estimate(Mdl,rets1(:,i));
        pqcoef(2*i-1,1:(1+pq(1)))=info.X(1:(1+pq(1)));
        for kk=1:(1+pq(1))
            pqcoef(2*i,kk)=sqrt(EstParamCov(kk,kk));
        end
        pqcoef(2*i-1,8:10)=info.X(end-2:end);
        for kk=1:3
            pqcoef(2*i,7+kk)=sqrt(EstParamCov(length(EstParamCov)-3+kk,length(EstParamCov)-3+kk));
        end
        pqcoef((2*i-1):(2*i),11)= outSKEWTandstd((2*i-1):(2*i),1);
        pqcoef((2*i-1):(2*i),12)= outSKEWTandstd((2*i-1):(2*i),2);
        pqcoef(2*i-1,13)=-logL;
        
    elseif pq(1)==0 & pq(2)>0
        Mdl = arima('MALags',[1:pq(2)],'Variance',garch(1,1));
        [EstMdl,EstParamCov,logL,info]  = estimate(Mdl,rets1(:,i));
        pqcoef(2*i-1,1)=info.X(1);
        pqcoef(2*i,1)=sqrt(EstParamCov(1,1));
        pqcoef(2*i-1,5:(4+pq(2)))=info.X(2:(pq(2)+1));
        for kk=2:(1+pq(2))
            pqcoef(2*i,3+kk)=sqrt(EstParamCov(kk,kk));
        end
        pqcoef(2*i-1,8:10)=info.X(end-2:end);
        for kk=1:3
            pqcoef(2*i,7+kk)=sqrt(EstParamCov(length(EstParamCov)-3+kk,length(EstParamCov)-3+kk));
        end
        pqcoef((2*i-1):(2*i),11)= outSKEWTandstd((2*i-1):(2*i),1);
        pqcoef((2*i-1):(2*i),12)= outSKEWTandstd((2*i-1):(2*i),2);
        pqcoef(2*i-1,13)=-logL;
        
        
    elseif pq(1)>0 & pq(2)>0
        Mdl = arima('ARLags',[1:pq(1)],'MALags',[1:pq(2)],'Variance',garch(1,1));
        [EstMdl,EstParamCov,logL,info]  = estimate(Mdl,rets1(:,i));
        
        pqcoef(2*i-1,1:(1+pq(1)))=info.X(1:(1+pq(1)));
        for kk=1:(1+pq(1))
            pqcoef(2*i,kk)=sqrt(EstParamCov(kk,kk));
        end
        pqcoef(2*i-1,5:(4+pq(2)))=info.X((2+pq(1)):(pq(2)+1+pq(1)));
        for kk=2:(1+pq(2))
            pqcoef(2*i,3+kk)=sqrt(EstParamCov(pq(1)+kk,pq(1)+kk));
        end
        pqcoef(2*i-1,8:10)=info.X(end-2:end);
        for kk=1:3
            pqcoef(2*i,7+kk)=sqrt(EstParamCov(length(EstParamCov)-3+kk,length(EstParamCov)-3+kk));
        end
        pqcoef((2*i-1):(2*i),11)= outSKEWTandstd((2*i-1):(2*i),1);
        pqcoef((2*i-1):(2*i),12)= outSKEWTandstd((2*i-1):(2*i),2);
        pqcoef(2*i-1,13)=-logL;
 
    end
    [h,pValue] = lbqtest(stdresids(:,i),'lags',10);
    pqcoef(2*i-1,14)=pValue;     %%p-value of Ljung–Box test
    [h,pValue] = lbqtest(stdresids(:,i).^2,'lags',10);
    pqcoef(2*i-1,15)=pValue;     %%%the Ljung–Box test of serial autocorrelation of squared errors
    [h,pValue,stat] = archtest(stdresids(:,i),'Lags',10,'Alpha',0.01);
    pqcoef(2*i-1,16)=pValue;     %%%ARCHtest 
    if i==1 | i==17 | i==18 |i==19   %%%这几个序列出现了停牌现象
        numb1=find(rets1(:,i)~=0);
        [h,pValue,stat] = kstest(stdresids(numb1,i),[stdresids(numb1,i)  skewtdis_cdf(stdresids(numb1,i),outSKEWT(i,1)+1,outSKEWT(i,2))],0.01);
        pqcoef(2*i-1,17)=pValue;
    else 
        [h,pValue,stat] = kstest(stdresids(:,i),[stdresids(:,i)  skewtdis_cdf(stdresids(:,i),outSKEWT(i,1),outSKEWT(i,2))],0.01);
        pqcoef(2*i-1,17)=pValue; 
    end
    
    [i toc]
end


pqcoefsta=cell(size(pqcoef,1),size(pqcoef,2));
for i=1:(size(pqcoef,1)/2)
    for j=1:(size(pqcoef,2)-5)
        pa=pqcoef(2*i-1,j);
        pb=pqcoef(2*i,j);
       
        if ~isnan(pa)
             pstd=abs(pa/pb);
            if pstd>=norminv(0.995)
                pstd1=strcat(num2str(pa,'%.3f'),'***');
                pstd2=strcat('(',num2str(pb,'%.3f'),')');
                pqcoefsta(2*i-1,j)=cellstr(pstd1);
                pqcoefsta(2*i,j)=cellstr(pstd2);
            elseif pstd<norminv(0.995) &pstd>=norminv(0.975)
                pstd1=strcat(num2str(pa,'%.3f'),'**');
                pstd2=strcat('(',num2str(pb,'%.3f'),')');
                pqcoefsta(2*i-1,j)=cellstr(pstd1);
                pqcoefsta(2*i,j)=cellstr(pstd2);
            elseif pstd<norminv(0.975)&pstd>=norminv(0.95)
                pstd1=strcat(num2str(pa,'%.3f'),'*');
                pstd2=strcat('(',num2str(pb,'%.3f'),')');
                pqcoefsta(2*i-1,j)=cellstr(pstd1);
                pqcoefsta(2*i,j)=cellstr(pstd2);
            elseif pstd<norminv(0.95)
                pstd1=strcat(num2str(pa,'%.3f'));
                pstd2=strcat('(',num2str(pb,'%.3f'),')');
                pqcoefsta(2*i-1,j)=cellstr(pstd1);
                pqcoefsta(2*i,j)=cellstr(pstd2);
            end
        end
    end
    for k=13:17
        pstd1=strcat(num2str(pqcoef(2*i-1,k),'%.3f'));
        pqcoefsta(2*i-1,k)=cellstr(pstd1);
    end
      
end


info.rnames = strvcat('石油期货','中国石油');


info.cnames = strvcat('.','phi0','phi1','phi2','phi3','theta1','theta2','theta3','w','alpha','beta','υ','λ','logL','LB','LB2','ARCH','K-S');
pqcoefsta2=cell(size(rets1,2)*2+1,size(pqcoef,2)+1);
pqcoefsta2(1,:)=cellstr(info.cnames)';
for i=2:(size(rets1,2)*2+1)
    for j=2:(size(pqcoef,2)+1)
        pqcoefsta2(i,j)=pqcoefsta(i-1,j-1);
    end
end

for i=1:size(rets1,2)
    pqcoefsta2(2*i,1)=cellstr(info.rnames(i,:));
end

%%%导出arma-garch结果

xlswrite('arma-garch结果1.xlsx',pqcoefsta2)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATING SOME CONSTANT COPULAS


options = optimset('Display','off','TolCon',10^-12,'TolFun',10^-4,'TolX',10^-6,'DiffMaxChange',Inf,'DiffMinChange',0,'Algorithm','active-set');

thetaALL=nan(2,5,7,3,size(rets1,2),size(rets1,2)); %2：第一行为系数值，第二行为系数标准差；5：前两列为对应系数，第三列为logL，第四列为AIC，第五列为BIC
                                                   %7：为[copula model]的类型；3：对应分布类型，[EDF or SKEWT marginal dist or T marginal dist]

LLALL = nan(7,3,size(rets1,2),size(rets1,2));   % logL; [ copula model] ; [EDF or SKEWT marginal dist or T dist]
tauLUall = nan(7,2,3,size(rets1,2),size(rets1,2));   % tail dependence
opt_copula=nan(size(rets1,2),size(rets1,2),3,2);  %%%%optimal copulas; 3:对应三种分布；2,：分别根据aic和bic信息准则判断
opt_aicbic=nan(size(rets1,2),size(rets1,2),3,2);  %%%%optimal copulas对应的最小的aic和bic

tic;
for i=1:(size(rets1,2)-1)
    for j=(i+1):size(rets1,2)
        for uu=2:3;  %  uu=2 uses SKEWT marginal dist,uu=3 uses T marginal dist
            u = Uall(:,i,uu);
            v = Uall(:,j,uu);
            % 1. Normal Copula
            kappa1 = corrcoef12(norminv(u),norminv(v));

            scoresc = LLgrad_1('NormalCopula_CLa',kappa1,[u v]);
            Hc1 = hessian('NormalCopula_CL',kappa1,[u v])/T;  % used in naive VCV matrix for copula params
            Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa1, sqrt(Vc1naive2)]; 
            LL1 = NormalCopula_CL(kappa1,[u,v]);
            [aic bic]=aicbic(-LL1,1,length(rets1));
            thetaALL(:,1,1,uu,i,j) =[kappa1, sqrt(Vc1naive2)]';
            thetaALL(:,3,1,uu,i,j) =-LL1;
            thetaALL(:,4,1,uu,i,j) =aic;
            thetaALL(:,5,1,uu,i,j) =bic; 

            % 2. Clayton's copula
            lower = 0.0001;
            theta0 = 1;
            [ kappa2 LL2] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);
            scoresc = LLgrad_1('claytonCLa',kappa2,[u v]);
            Hc1 = hessian('claytonCL',kappa2,[u v])/T;  % used in naive VCV matrix for copula params
            Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa2, sqrt(Vc1naive2)]; 
            [aic bic]=aicbic(-LL2,1,length(rets1));
            thetaALL(:,1,2,uu,i,j) =[kappa2, sqrt(Vc1naive2)]';
            thetaALL(:,3,2,uu,i,j) =-LL2;
            thetaALL(:,4,2,uu,i,j) =aic;
            thetaALL(:,5,2,uu,i,j) =bic; 

            % 3. Rotated Clayton copula (with tail dep in upper tail instead of lower)
            lower = 0.0001;
            theta0 = 1;
           [ kappa3 LL3] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,1-[u,v]);
            scoresc = LLgrad_1('claytonCLa',kappa3,1-[u v]);
            Hc1 = hessian('claytonCL',kappa3,1-[u v])/T;  % used in naive VCV matrix for copula params
            Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa3, sqrt(Vc1naive2)]; 
            [aic bic]=aicbic(-LL3,1,length(rets1));
            thetaALL(:,1,3,uu,i,j) =[kappa3, sqrt(Vc1naive2)]';
            thetaALL(:,3,3,uu,i,j) =-LL3;
            thetaALL(:,4,3,uu,i,j) =aic;
            thetaALL(:,5,3,uu,i,j) =bic; 

           
         
            % 4. Gumbel copula
            lower = 1.1;
            upper = 5;
            theta0 = 2;
            [ kappa4 LL4] = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
            scoresc = LLgrad_1('gumbelCLa',kappa4,[u v]);
            Hc1 = hessian('gumbelCL',kappa4,[u v])/T;  % used in naive VCV matrix for copula params
            Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa4, sqrt(Vc1naive2)]; 
            [aic bic]=aicbic(-LL4,1,length(rets1));
            thetaALL(:,1,4,uu,i,j) =[kappa4, sqrt(Vc1naive2)]';
            thetaALL(:,3,4,uu,i,j) =-LL4;
            thetaALL(:,4,4,uu,i,j) =aic;
            thetaALL(:,5,4,uu,i,j) =bic; 

             % 5. Rotated Gumbel copula
            lower = 1.1;
            upper = 5;
            theta0 = 2;
            [ kappa5 LL5] = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,1-[u,v]);
            scoresc = LLgrad_1('gumbelCLa',kappa5,1-[u v]);
            Hc1 = hessian('gumbelCL',kappa5,1-[u v])/T;  % used in naive VCV matrix for copula params
            Vc1naive2 = --inv(Hc1)/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa5, sqrt(Vc1naive2)]; 
            [aic bic]=aicbic(-LL5,1,length(rets1));
            thetaALL(:,1,5,uu,i,j) =[kappa5, sqrt(Vc1naive2)]';
            thetaALL(:,3,5,uu,i,j) =-LL5;
            thetaALL(:,4,5,uu,i,j) =aic;
            thetaALL(:,5,5,uu,i,j) =bic; 
            
            
    
            % 6. Symmetrised Joe-Clayton copula
           lower = [0 , 0 ];
           upper = [ 1 , 1];
           theta0 = [0.4;0.4];
           [ kappa6 LL6] = fmincon('sym_jc_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
            scoresc = LLgrad_1('sym_jc_CLa',kappa6,[u v]);
            Hc1 = hessian('sym_jc_CL',kappa6,[u v])/T;  % used in naive VCV matrix for copula params
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;
             % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
            if Vc2(1,1)<0
               i
               j
            end


            [kappa6, sqrt(diag(Vc2))]; 
            [aic bic]=aicbic(-LL6,2,length(rets1));   %两个参数
            
            thetaALL(:,1:2,6,uu,i,j) =[kappa6, sqrt(diag(Vc2))]';
            thetaALL(:,3,6,uu,i,j) =-LL6;
            thetaALL(:,4,6,uu,i,j) =aic;
            thetaALL(:,5,6,uu,i,j) =bic;  
           
           
           % 7. Student's t copula  (estimating nu_inv rather than nu)
           lower = [-0.9 , 1/100 ];
           upper = [ 0.9 , 1/2.1 ];
           theta0 = [kappa1;10];
           [kappa7 LL7] = fmincon('tcopulaCL2',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
            scoresc = LLgrad_1('tcopulaCL2a',kappa7,[u v]);
            Hc1 = hessian('tcopulaCL2',kappa7,[u v])/T;  % used in naive VCV matrix for copula params
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive VCV matrix for copula params. ignoring estimation error from the marginal distributions
      
           if Vc2(1,1)<0 | Vc2(2,2)<0
                  stu71=i
                  stu72=j
            end
            [kappa7, sqrt(diag(Vc2))]; 
            [aic bic]=aicbic(-LL7,2,length(rets1));   %两个参数
            thetaALL(:,1:2,7,uu,i,j) =[kappa7, sqrt(diag(Vc2))]';
            thetaALL(:,3,7,uu,i,j) =-LL7;
            thetaALL(:,4,7,uu,i,j) =aic;
            thetaALL(:,5,7,uu,i,j) =bic;    
            if LL6==-inf
                LL6=max([LL1;LL2;LL3;LL4;LL5;LL6;LL7]);
            end
            

           LLALL(:,uu,i,j) = -[LL1;LL2;LL3;LL4;LL5;LL6;LL7];
           copula_params = [ones(5,1);2;2];
           [aic bic]=aicbic(LLALL(:,uu,i,j),copula_params,length(rets1));
           opt_copula(i,j,uu,1) = find(aic==min(aic));
           opt_aicbic(i,j,uu,1) = min(aic);
           opt_copula(i,j,uu,2) = find(bic==min(bic));
           opt_aicbic(i,j,uu,2) = min(bic);
           
           
            % tail dependence implied by each of these copulas
            tauLUall(1,:,uu,i,j) = [0,0];                 % Normal copula has zero tail dependence
            tauLUall(2,:,uu,i,j) = [2^(-1/kappa2),0];     % Clayton copula has zero upper tail dependence
            tauLUall(3,:,uu,i,j) = [0,2^(-1/kappa3)];     % Rotated Clayton copula has zero lower tail dependence
            tauLUall(4,:,uu,i,j) = [0,2-2^(1/kappa4)];    % Gumbel copula has zero lower tail dependence
            tauLUall(5,:,uu,i,j) = [2-2^(1/kappa5),0];    % Rotated Gumbel copula has zero upper tail dependence
            
            tauLUall(6,:,uu,i,j) = kappa6([2,1])';               % SJC copula parameters are the tail dependence coefficients, but in reverse order.
            tauLUall(7,:,uu,i,j) = ones(1,2)*2*tdis_cdf(-sqrt((kappa7(2)+1)*(1-kappa7(1))/(1+kappa7(1))),kappa7(2)+1);  % Student's t copula has symmetric tail dependence
        end
       
    end
    [i,toc]
end

 
xlswrite('静态最优Copula2.xlsx',opt_copula(:,:,2,1))   %%%2表示偏t分布，根据AIC准则的最优Copula
length(find(opt_copula(:,:,2,1)==7));
length(find(opt_copula(:,:,2,1)==6));
length(find(opt_copula(:,:,2,1)==5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ESTIMATING SOME dynamic COPULAS,五种动态copula

dynamic_thetaALL =nan(2,9,7,3,size(rets1,2),size(rets1,2)); %2：第一行为系数值，第二行为系数标准差；9：前六列为对应系数，第七列为logL，第八列为AIC，第九列为BIC
                                                   %7：为[copula model]的类型；3：对应分布类型，[EDF or SKEWT marginal dist or T marginal dist]
dynamic_LLALL = nan(7,3,size(rets1,2),size(rets1,2));   %  % logL; [ copula model] ; [EDF or SKEWT marginal dist or T dist]

dynamic_rho=nan(length(rets1),2,7,3,size(rets1,2),size(rets1,2));   % 接收动态的rho值,2表示有可能存在两个rho; 7:表示动态copula模型个数；

dynamic_opt_copula=nan(size(rets1,2),size(rets1,2),3,2);  %%%%optimal copulas; 3:对应三种分布；2,：分别根据aic和bic信息准则判断
dynamic_opt_aicbic=nan(size(rets1,2),size(rets1,2),3,2);  %%%%optimal copulas对应的最小的aic和bic

tic;
for i=1:(size(rets1,2)-1)
    for j=(i+1):size(rets1,2)
        for uu=2:3;  %  uu=2 uses SKEWT marginal dist,uu=3 uses T marginal dist
            u = Uall(:,i,uu);
            v = Uall(:,j,uu);
            %%%%%%11 Time-varying normal Copula
            kappa1 = corrcoef12(norminv(u),norminv(v));
                
            lower = -5*ones(3,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
            upper = 5*ones(3,1);
            theta0 = [log((1+kappa1)/(1-kappa1));0;0];
            [kappa11 LL11] = fmincon('bivnorm_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa1);
            [LL11, rho11] = bivnorm_tvp1_CL(kappa11,[u,v],kappa1);
            
            scoresc = LLgrad_1('bivnorm_tvp1_CLa',kappa11,[u,v],kappa1);
            Hc1 =     hessian('bivnorm_tvp1_CL', kappa11,[u,v],kappa1)/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa11, sqrt(diag(Vc2))];
            [aic bic]=aicbic(-LL11,3,length(rets1));
            dynamic_thetaALL(:,1:3,1,uu,i,j) =[kappa11, sqrt(diag(Vc2))]';
            dynamic_thetaALL(:,7,1,uu,i,j)   = -LL11;
            dynamic_thetaALL(:,8,1,uu,i,j)   = aic;
            dynamic_thetaALL(:,9,1,uu,i,j)   = bic;
            dynamic_rho(:,1,1,uu,i,j)=rho11;
  
            % 12. Time-varying clayton copula
            lower = 0.0001;
            theta0 = 1;
            [ kappa2 LL2] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,[u,v]);
            
            lower = -5*ones(3,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
            upper =  5*ones(3,1);
            theta0 = [sqrt(kappa2);0;0];
            [ kappa12 LL12] = fmincon('clayton_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa2);
            [LL12, rho12] = clayton_tvp1_CL(kappa12,[u,v],kappa2);

            scoresc = LLgrad_1('clayton_tvp1_CLa',kappa12,[u,v],kappa2);
            Hc1 =     hessian('clayton_tvp1_CL', kappa12,[u,v],kappa2)/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa12, sqrt(diag(Vc2)) ];
            [aic bic]=aicbic(-LL12,3,length(rets1));
            dynamic_thetaALL(:,1:3,2,uu,i,j) =[kappa12, sqrt(diag(Vc2))]';
            dynamic_thetaALL(:,7,2,uu,i,j)   = -LL12;
            dynamic_thetaALL(:,8,2,uu,i,j)   = aic;
            dynamic_thetaALL(:,9,2,uu,i,j)   = bic;
            dynamic_rho(:,1,2,uu,i,j)=rho12;
            
            % 13. Time-varying Rotated clayton copula
            lower = 0.0001;
            theta0 = 1;
            [ kappa3 LL3] = fmincon('claytonCL',theta0,[],[],[],[],lower,[],[],options,1-[u,v]);
            lower = -5*ones(3,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
            upper =  5*ones(3,1);
            theta0 = [sqrt(kappa3);0;0];
            [kappa13 LL13] = fmincon('clayton_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,1-[u,v],kappa3);
            [LL13, rho13] = clayton_tvp1_CL(kappa13,1-[u,v],kappa3);

            scoresc = LLgrad_1('clayton_tvp1_CLa',kappa13,1-[u,v],kappa3);
            Hc1 =     hessian('clayton_tvp1_CL', kappa13,1-[u,v],kappa3)/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa13, sqrt(diag(Vc2)) ];
            [aic bic]=aicbic(-LL13,3,length(rets1));
            dynamic_thetaALL(:,1:3,3,uu,i,j) =[kappa13, sqrt(diag(Vc2))]';
            dynamic_thetaALL(:,7,3,uu,i,j)   = -LL13;
            dynamic_thetaALL(:,8,3,uu,i,j)   = aic;
            dynamic_thetaALL(:,9,3,uu,i,j)   = bic;
            dynamic_rho(:,1,3,uu,i,j)=rho13; 
            
            % 14. Time-varying  Gumbel copula  
            lower = 1.1;
            upper = 5;
            theta0 = 2;
            [ kappa4 LL4] = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,[u,v]);
            lower = -5*ones(3,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
            upper =  5*ones(3,1);
            theta0 = [sqrt(kappa4-1);0;0];
            [kappa14 LL14] = fmincon('Gumbel_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa4);
            [LL14, rho14] = Gumbel_tvp1_CL(kappa14,[u,v],kappa4);

            scoresc = LLgrad_1('Gumbel_tvp1_CLa',kappa14,[u,v],kappa4);
            Hc1 =     hessian('Gumbel_tvp1_CL', kappa14,[u,v],kappa4)/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa14, sqrt(diag(Vc2))];
            [aic bic]=aicbic(-LL14,3,length(rets1));
            dynamic_thetaALL(:,1:3,4,uu,i,j) =[kappa14, sqrt(diag(Vc2))]';
            dynamic_thetaALL(:,7,4,uu,i,j)   = -LL14;
            dynamic_thetaALL(:,8,4,uu,i,j)   = aic;
            dynamic_thetaALL(:,9,4,uu,i,j)   = bic;
            dynamic_rho(:,1,4,uu,i,j)=rho14; 
            
            % 15. Time-varying Rotated Gumbel copula
            lower = 1.1;
            upper = 5;
            theta0 = 2;
            [kappa5 LL5] = fmincon('gumbelCL',theta0,[],[],[],[],lower,upper,[],options,1-[u,v]);          
            
            lower = -5*ones(3,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
            upper =  5*ones(3,1);
            theta0 = [sqrt(kappa5-1);0;0];
            [kappa15 LL15] = fmincon('Gumbel_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,1-[u,v],kappa5);
            [LL15, rho15] = Gumbel_tvp1_CL(kappa15,1-[u,v],kappa5);

            scoresc = LLgrad_1('Gumbel_tvp1_CLa',kappa15,1-[u,v],kappa5);
            Hc1 =     hessian('Gumbel_tvp1_CL', kappa15,1-[u,v],kappa5)/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa15, sqrt(diag(Vc2)) ];
            
            [aic bic]=aicbic(-LL15,3,length(rets1));
            dynamic_thetaALL(:,1:3,5,uu,i,j) =[kappa15, sqrt(diag(Vc2))]';
            dynamic_thetaALL(:,7,5,uu,i,j)   = -LL15;
            dynamic_thetaALL(:,8,5,uu,i,j)   = aic;
            dynamic_thetaALL(:,9,5,uu,i,j)   = bic;
            dynamic_rho(:,1,5,uu,i,j)=rho15;    
            
            % 16. Time-varying SJC copula
            lower = [0 , 0 ];
            upper = [ 1 , 1];
            theta0 = [0.4;0.4];
            [ kappa6 LL6] = fmincon('sym_jc_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v]); 
            lower = -25*ones(6,1);  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
            upper =  25*ones(6,1);
            theta0 = [log(kappa6(1)/(1-kappa6(1)));0;0;log(kappa6(2)/(1-kappa6(2)));0;0];
            [ kappa16 LL16] = fmincon('sym_jc_tvp_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v],kappa6);
            [ LL16 tauU16 tauL16] = sym_jc_tvp_CL(kappa16,[u,v],kappa6);

            scoresc = LLgrad_1('sym_jc_tvp_CLa',kappa16,[u,v],kappa6);
            Hc1 =     hessian('sym_jc_tvp_CL', kappa16,[u,v],kappa6)/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
            Bc = newey_west(scoresc,floor(T^(1/3)));
            
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
            
            [kappa16, sqrt(diag(Vc2)) ];

            [aic bic]=aicbic(-LL16,3,length(rets1));
            dynamic_thetaALL(:,1:6,6,uu,i,j) =[kappa16, sqrt(diag(Vc2))]';
            dynamic_thetaALL(:,7,6,uu,i,j)   = -LL16;
            dynamic_thetaALL(:,8,6,uu,i,j)   = aic;
            dynamic_thetaALL(:,9,6,uu,i,j)   = bic;
            dynamic_rho(:,1,6,uu,i,j)=tauU16;      
            dynamic_rho(:,2,6,uu,i,j)=tauL16;

            % 17 Time-varying students't copula
            lower = [-0.9 , 1/100 ];
            upper = [ 0.9 , 1/2.1 ];
            theta0 = [kappa1;10];
            [kappa7 LL7] = fmincon('tcopulaCL2',theta0,[],[],[],[],lower,upper,[],options,[u,v]);

            lower = [-5*ones(3,1);1];  % in theory there are no constraints, but setting loose constraints sometimes helps in the numerical optimisation
            upper = [5*ones(3,1);100];
            theta0 = [log( (0.9999+kappa7(1))/(0.9999-kappa7(1)) )*(1-0.05-0.8);0.05;0.8;1/kappa7(2)];
            [ kappa17 LL17] = fmincon('tcopula_tvp1_CL',theta0,[],[],[],[],lower,upper,[],options,[u,v],[kappa7(1);1/kappa7(2)] );
            [LL17, rho17] = tcopula_tvp1_CL(kappa17,[u,v],[kappa7(1);1/kappa7(2)]);

            scoresc = LLgrad_1('tcopula_tvp1_CLa',kappa17,[u,v],[kappa7(1);1/kappa7(2)]);
            Hc1 =     hessian('tcopula_tvp1_CL', kappa17,[u,v],[kappa7(1);1/kappa7(2)])/T;  % used in naive VCV matrix for copula params [h = -eps.^(1/3)*max(abs(x),1e-4);]
            Bc = newey_west(scoresc,floor(T^(1/3)));
            Vc2 = inv(Hc1)*Bc*(inv(Hc1)')/T;  % naive sandwich VCV matrix for copula params. ignoring estimation error from the marginal distributions
            [kappa17, sqrt(diag(Vc2)) ];
 
            [aic bic]=aicbic(-LL17,3,length(rets1));
            dynamic_thetaALL(:,1:4,7,uu,i,j) =[kappa17, sqrt(diag(Vc2))]';
            dynamic_thetaALL(:,7,7,uu,i,j)   = -LL17;
            dynamic_thetaALL(:,8,7,uu,i,j)   = aic;
            dynamic_thetaALL(:,9,7,uu,i,j)   = bic;
            dynamic_rho(:,1,7,uu,i,j)=rho17;


            if LL16==-inf
                LL16=max([LL11;LL12;LL13;LL14;LL15;LL16;LL17]);
            end

            dynamic_LLALL(:,uu,i,j) = -[LL11;LL12;LL13;LL14;LL15;LL16;LL17];
            copula_params = [3;3;3;3;3;6;4];
            [aic bic]=aicbic(dynamic_LLALL(:,uu,i,j),copula_params,length(rets1));
            dynamic_opt_copula(i,j,uu,1) = find(aic==min(aic));
            dynamic_opt_aicbic(i,j,uu,1) = min(aic);
            dynamic_opt_copula(i,j,uu,2) = find(bic==min(bic));
            dynamic_opt_aicbic(i,j,uu,2) = min(bic);

            toc
        end
        [i j toc]  
    end
end

xlswrite('动态最优copula2.xlsx',dynamic_opt_copula(:,:,2,1))
length(find(dynamic_opt_copula(:,:,2,1)==7));
length(find(dynamic_opt_copula(:,:,2,1)==6));
length(find(dynamic_opt_copula(:,:,2,1)==5));

%%%%%%极大似然估计值
a1=dynamic_LLALL(:,2,:,:);
b1=permute(a1,[1 3 4 2]);

[c1 c2 c3]=find(b1==inf);

[i,j,z]=ind2sub(size(b1),find(b1==inf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%VALUE-AT-RISK （VaR）and conditional VaR (CoVaR) and the delta CoVaR FROM COPULA-BASED MODELS
 
QQ = [0.01;0.05;0.5;0.95;0.99];   %%%%显著性水平


muhat = rets1 - resids;  % conditional mean

outVAR = nan(T,size(rets1,2),length(QQ),3);  %%%对应downside和Upside两个VaR,3:对应对应分布类型，[EDF or SKEWT marginal dist or T marginal dist]

tic;
for ii=1:(size(rets1,2))
    for tt=1:length(rets1)
        uandlamdai=outSKEWT(ii,:);  %%%第i列偏t分布对应的参数自由度nu,lambda
        uofstudentst=tFreedom(ii);  %%%第i列t分布对应的参数自由度nu
        outVAR(tt,ii,:,2) =muhat(tt,ii) + sqrt(hhat_opt(tt,ii))*skewtdis_inv(QQ, uandlamdai(1), uandlamdai(2));
        outVAR(tt,ii,:,3) =muhat(tt,ii) + sqrt(hhat_opt(tt,ii))*tinv(QQ,uofstudentst);
    end
    [ii toc]
end
toc


%%%%%%
figure(13);
plot(datet, outVAR(:,1,2,3),'LineWidth',2)  %3对应t分布
hold on
plot(datet, outVAR(:,1,4,3),'LineWidth',2)


%%%%%%%%根据copula函数和CoVaR和DeltaCoVaR
outCoVAR=nan(T,size(rets1,2),size(rets1,2),length(QQ),3,7); %%%3:对应对应分布类型，[EDF or SKEWT marginal dist or T marginal dist];
                                                            %%%7:表示7中时变copula模型
 
outDeltaCoVAR=nan(T,size(rets1,2),size(rets1,2),length(QQ),3,7); %%%3:对应对应分布类型，[EDF or SKEWT marginal dist or T marginal dist];
                                                            %%%7:表示7中时变copula模型
 
tic
for i=1:(size(rets1,2)-1)
    for j=(i+1):size(rets1,2)
        for tt=1:length(rets1)
            for qq=1:length(QQ)
                %%%uu=2
                uu=2;
                vari=outVAR(tt,i,qq,uu);
                varj=outVAR(tt,j,qq,uu);
                uandlamdai=outSKEWT(i,:);  %%%第i列偏t分布对应的参数自由度nu,lambda
                uandlamdaj=outSKEWT(j,:);  %%%第j列偏t分布对应的参数nu,lambda
                if QQ(qq)<0.5
                    %%%%%%11 Time-varying normal Copula
                    rho11=dynamic_rho(tt,1,1,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u QQ(qq)],rho11)-QQ(qq)^2,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,1)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,1)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u 0.5],rho11)-0.5*QQ(qq),[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,1)=(outCoVAR(tt,i,j,qq,uu,1)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,1)=(outCoVAR(tt,j,i,qq,uu,1)-covarmedxj)/covarmedxj;    

                    % 12. Time-varying clayton copula
                    rho12=dynamic_rho(tt,1,2,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho12)-QQ(qq)^2,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,2)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,2)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho12)-0.5*QQ(qq),[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,2)=(outCoVAR(tt,i,j,qq,uu,2)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,2)=(outCoVAR(tt,j,i,qq,uu,2)-covarmedxj)/covarmedxj;    

                     % 13. Time-varying Rotated clayton copula
                    rho13=dynamic_rho(tt,1,3,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho13)-QQ(qq)^2,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,3)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,3)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho13)-0.5*QQ(qq),[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,3)=(outCoVAR(tt,i,j,qq,uu,3)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,3)=(outCoVAR(tt,j,i,qq,uu,3)-covarmedxj)/covarmedxj;    

                     % 14. Time-varying  Gumbel copula
                    rho14=dynamic_rho(tt,1,4,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho14)-QQ(qq)^2,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,4)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,4)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho14)-0.5*QQ(qq),[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,4)=(outCoVAR(tt,i,j,qq,uu,4)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,4)=(outCoVAR(tt,j,i,qq,uu,4)-covarmedxj)/covarmedxj;    
                    
                     % 15. Time-varying Rotated Gumbel copula
                     rho15=dynamic_rho(tt,1,5,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho15)-QQ(qq)^2,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,5)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,5)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho15)-0.5*QQ(qq),[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,5)=(outCoVAR(tt,i,j,qq,uu,5)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,5)=(outCoVAR(tt,j,i,qq,uu,5)-covarmedxj)/covarmedxj;    
                    
                    % 16. Time-varying SJC copula
                    tauU12=dynamic_rho(tt,1,6,uu,i,j);  %%%Time-varying SJC copula
                    tauL12=dynamic_rho(tt,2,6,uu,i,j);
                    try
                       [sjc_u,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,QQ(qq),tauU12,tauL12)-QQ(qq)^2,[0 1]);
                       xi=skewtdis_inv(sjc_u, uandlamdai(1), uandlamdai(2));
                       xj=skewtdis_inv(sjc_u, uandlamdaj(1), uandlamdaj(2));
                       outCoVAR(tt,i,j,qq,uu,6)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;
                       outCoVAR(tt,j,i,qq,uu,6)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;
                    catch
                       outReporterrors(tt,i,j,qq,uu,1)=1;
                    end
                    [sjc_medu,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,0.5,tauU12,tauL12)-0.5*QQ(qq),[0 1]);    %%%VaR的置信水平为0.5时求得的CoVaR
                    medxi=skewtdis_inv(sjc_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(sjc_medu, uandlamdaj(1), uandlamdaj(2));
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                        
                    outDeltaCoVAR(tt,i,j,qq,uu,6)=(outCoVAR(tt,i,j,qq,uu,6)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,6)=(outCoVAR(tt,j,i,qq,uu,6)-covarmedxj)/covarmedxj;
                    
                     % 17 Time-varying students't copula
                     rho17=dynamic_rho(tt,1,7,uu,i,j);
                     theta8tv2= dynamic_thetaALL(1,4,7,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('t',[u QQ(qq)],rho17,theta8tv2)-QQ(qq)^2,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,7)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,7)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('t',[u 0.5],rho17,theta8tv2)-0.5*QQ(qq),[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,7)=(outCoVAR(tt,i,j,qq,uu,7)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,7)=(outCoVAR(tt,j,i,qq,uu,7)-covarmedxj)/covarmedxj;    
                    
                 else
                     alpha0=1-QQ(qq);
                     beta0=1-QQ(qq);
                     %%%%%%11 Time-varying normal Copula
                     rho11=dynamic_rho(tt,1,1,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u QQ(qq)],rho11)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,1)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,1)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u 0.5],rho11)-u+0.5-0.5*beta0,[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,1)=(outCoVAR(tt,i,j,qq,uu,1)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,1)=(outCoVAR(tt,j,i,qq,uu,1)-covarmedxj)/covarmedxj;    

                    % 12. Time-varying clayton copula
                    rho12=dynamic_rho(tt,1,2,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho12)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,2)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,2)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho12)-u+0.5-0.5*beta0,[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,2)=(outCoVAR(tt,i,j,qq,uu,2)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,2)=(outCoVAR(tt,j,i,qq,uu,2)-covarmedxj)/covarmedxj;    

                     % 13. Time-varying Rotated clayton copula
                    rho13=dynamic_rho(tt,1,3,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho13)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,3)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,3)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho13)-u+0.5-0.5*beta0,[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,3)=(outCoVAR(tt,i,j,qq,uu,3)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,3)=(outCoVAR(tt,j,i,qq,uu,3)-covarmedxj)/covarmedxj;    

                     % 14. Time-varying  Gumbel copula
                    rho14=dynamic_rho(tt,1,4,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho14)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,4)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,4)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho14)-u+0.5-0.5*beta0,[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,4)=(outCoVAR(tt,i,j,qq,uu,4)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,4)=(outCoVAR(tt,j,i,qq,uu,4)-covarmedxj)/covarmedxj;    
                    
                     % 15. Time-varying Rotated Gumbel copula
                     rho15=dynamic_rho(tt,1,5,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho15)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,5)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,5)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho15)-u+0.5-0.5*beta0,[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,5)=(outCoVAR(tt,i,j,qq,uu,5)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,5)=(outCoVAR(tt,j,i,qq,uu,5)-covarmedxj)/covarmedxj;    
                    
                    % 16. Time-varying SJC copula
                    tauU12=dynamic_rho(tt,1,6,uu,i,j);  %%%Time-varying SJC copula
                    tauL12=dynamic_rho(tt,2,6,uu,i,j);
                    try
                       [sjc_u,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,QQ(qq),tauU12,tauL12)-u+alpha0-alpha0*beta0,[0 1]);
                       xi=skewtdis_inv(sjc_u, uandlamdai(1), uandlamdai(2));
                       xj=skewtdis_inv(sjc_u, uandlamdaj(1), uandlamdaj(2));
                       outCoVAR(tt,i,j,qq,uu,6)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;
                       outCoVAR(tt,j,i,qq,uu,6)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;
                    catch
                       outReporterrors(tt,i,j,qq,uu,1)=1;
                    end
                    [sjc_medu,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,0.5,tauU12,tauL12)-u+0.5-0.5*beta0,[0 1]);    %%%VaR的置信水平为0.5时求得的CoVaR
                    medxi=skewtdis_inv(sjc_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(sjc_medu, uandlamdaj(1), uandlamdaj(2));
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                        
                    outDeltaCoVAR(tt,i,j,qq,uu,6)=(outCoVAR(tt,i,j,qq,uu,6)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,6)=(outCoVAR(tt,j,i,qq,uu,6)-covarmedxj)/covarmedxj;
                    
                     % 17 Time-varying students't copula
                     rho17=dynamic_rho(tt,1,7,uu,i,j);
                     theta8tv2= dynamic_thetaALL(1,4,7,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('t',[u QQ(qq)],rho17,theta8tv2)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=skewtdis_inv(fu, uandlamdai(1), uandlamdai(2));
                    xj=skewtdis_inv(fu, uandlamdaj(1), uandlamdaj(2));
                    outCoVAR(tt,i,j,qq,uu,7)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,7)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('t',[u 0.5],rho17,theta8tv2)-u+0.5-0.5*beta0,[0 1]);
                    medxi=skewtdis_inv(f_medu, uandlamdai(1), uandlamdai(2));
                    medxj=skewtdis_inv(f_medu, uandlamdaj(1), uandlamdaj(2)); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,7)=(outCoVAR(tt,i,j,qq,uu,7)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,7)=(outCoVAR(tt,j,i,qq,uu,7)-covarmedxj)/covarmedxj;
                 end
                 

                 %%%uu=3
                 uu=3;
                 vari=outVAR(tt,i,qq,uu);
                 varj=outVAR(tt,j,qq,uu);
                 uandlamdai=tFreedom(i);  %%%第i列偏t分布对应的参数自由度nu,lambda
                 uandlamdaj=tFreedom(j);  %%%第j列偏t分布对应的参数nu,lambda
                 if QQ(qq)<0.5
                    %%%%%%11 Time-varying normal Copula
                    rho11=dynamic_rho(tt,1,1,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u QQ(qq)],rho11)-QQ(qq)^2,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,1)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,1)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u 0.5],rho11)-0.5*QQ(qq),[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,1)=(outCoVAR(tt,i,j,qq,uu,1)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,1)=(outCoVAR(tt,j,i,qq,uu,1)-covarmedxj)/covarmedxj;    

                    % 12. Time-varying clayton copula
                    rho12=dynamic_rho(tt,1,2,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho12)-QQ(qq)^2,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,2)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,2)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho12)-0.5*QQ(qq),[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,2)=(outCoVAR(tt,i,j,qq,uu,2)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,2)=(outCoVAR(tt,j,i,qq,uu,2)-covarmedxj)/covarmedxj;    

                     % 13. Time-varying Rotated clayton copula
                    rho13=dynamic_rho(tt,1,3,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho13)-QQ(qq)^2,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,3)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,3)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho13)-0.5*QQ(qq),[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,3)=(outCoVAR(tt,i,j,qq,uu,3)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,3)=(outCoVAR(tt,j,i,qq,uu,3)-covarmedxj)/covarmedxj;    

                     % 14. Time-varying  Gumbel copula
                    rho14=dynamic_rho(tt,1,4,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho14)-QQ(qq)^2,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,4)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,4)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho14)-0.5*QQ(qq),[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,4)=(outCoVAR(tt,i,j,qq,uu,4)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,4)=(outCoVAR(tt,j,i,qq,uu,4)-covarmedxj)/covarmedxj;    
                    
                     % 15. Time-varying Rotated Gumbel copula
                     rho15=dynamic_rho(tt,1,5,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho15)-QQ(qq)^2,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,5)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,5)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho15)-0.5*QQ(qq),[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,5)=(outCoVAR(tt,i,j,qq,uu,5)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,5)=(outCoVAR(tt,j,i,qq,uu,5)-covarmedxj)/covarmedxj;    
                    
                    % 16. Time-varying SJC copula
                    tauU12=dynamic_rho(tt,1,6,uu,i,j);  %%%Time-varying SJC copula
                    tauL12=dynamic_rho(tt,2,6,uu,i,j);
                    try
                       [sjc_u,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,QQ(qq),tauU12,tauL12)-QQ(qq)^2,[0 1]);
                       xi=tinv(sjc_u, uandlamdai);
                       xj=tinv(sjc_u, uandlamdaj);
                       outCoVAR(tt,i,j,qq,uu,6)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;
                       outCoVAR(tt,j,i,qq,uu,6)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;
                    catch
                       outReporterrors(tt,i,j,qq,uu,1)=1;
                    end
                    [sjc_medu,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,0.5,tauU12,tauL12)-0.5*QQ(qq),[0 1]);    %%%VaR的置信水平为0.5时求得的CoVaR
                    medxi=tinv(sjc_medu, uandlamdai);
                    medxj=tinv(sjc_medu, uandlamdaj);
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                        
                    outDeltaCoVAR(tt,i,j,qq,uu,6)=(outCoVAR(tt,i,j,qq,uu,6)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,6)=(outCoVAR(tt,j,i,qq,uu,6)-covarmedxj)/covarmedxj;
                    
                     % 17 Time-varying students't copula
                     rho17=dynamic_rho(tt,1,7,uu,i,j);
                     theta8tv2= dynamic_thetaALL(1,4,7,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('t',[u QQ(qq)],rho17,theta8tv2)-QQ(qq)^2,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,7)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,7)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('t',[u 0.5],rho17,theta8tv2)-0.5*QQ(qq),[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,7)=(outCoVAR(tt,i,j,qq,uu,7)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,7)=(outCoVAR(tt,j,i,qq,uu,7)-covarmedxj)/covarmedxj;    
                    
                 else
                     alpha0=1-QQ(qq);
                     beta0=1-QQ(qq);
                     %%%%%%11 Time-varying normal Copula
                     rho11=dynamic_rho(tt,1,1,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u QQ(qq)],rho11)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,1)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,1)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gaussian',[u 0.5],rho11)-u+0.5-0.5*beta0,[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,1)=(outCoVAR(tt,i,j,qq,uu,1)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,1)=(outCoVAR(tt,j,i,qq,uu,1)-covarmedxj)/covarmedxj;    

                    % 12. Time-varying clayton copula
                    rho12=dynamic_rho(tt,1,2,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho12)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,2)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,2)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho12)-u+0.5-0.5*beta0,[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,2)=(outCoVAR(tt,i,j,qq,uu,2)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,2)=(outCoVAR(tt,j,i,qq,uu,2)-covarmedxj)/covarmedxj;    

                     % 13. Time-varying Rotated clayton copula
                    rho13=dynamic_rho(tt,1,3,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u QQ(qq)],rho13)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,3)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,3)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Clayton',[u 0.5],rho13)-u+0.5-0.5*beta0,[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,3)=(outCoVAR(tt,i,j,qq,uu,3)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,3)=(outCoVAR(tt,j,i,qq,uu,3)-covarmedxj)/covarmedxj;    

                     % 14. Time-varying  Gumbel copula
                    rho14=dynamic_rho(tt,1,4,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho14)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,4)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,4)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho14)-u+0.5-0.5*beta0,[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,4)=(outCoVAR(tt,i,j,qq,uu,4)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,4)=(outCoVAR(tt,j,i,qq,uu,4)-covarmedxj)/covarmedxj;    
                    
                     % 15. Time-varying Rotated Gumbel copula
                     rho15=dynamic_rho(tt,1,5,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u QQ(qq)],rho15)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,5)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,5)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('Gumbel',[u 0.5],rho15)-u+0.5-0.5*beta0,[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,5)=(outCoVAR(tt,i,j,qq,uu,5)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,5)=(outCoVAR(tt,j,i,qq,uu,5)-covarmedxj)/covarmedxj;    
                    
                    % 16. Time-varying SJC copula
                    tauU12=dynamic_rho(tt,1,6,uu,i,j);  %%%Time-varying SJC copula
                    tauL12=dynamic_rho(tt,2,6,uu,i,j);
                    try
                       [sjc_u,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,QQ(qq),tauU12,tauL12)-u+alpha0-alpha0*beta0,[0 1]);
                       xi=tinv(sjc_u, uandlamdai);
                       xj=tinv(sjc_u, uandlamdaj);
                       outCoVAR(tt,i,j,qq,uu,6)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;
                       outCoVAR(tt,j,i,qq,uu,6)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;
                    catch
                       outReporterrors(tt,i,j,qq,uu,1)=1;
                    end
                    [sjc_medu,fval,exitflag]=fzero(@(u) sym_jc_cdf(u,0.5,tauU12,tauL12)-u+0.5-0.5*beta0,[0 1]);    %%%VaR的置信水平为0.5时求得的CoVaR
                    medxi=tinv(sjc_medu, uandlamdai);
                    medxj=tinv(sjc_medu, uandlamdaj);
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                        
                    outDeltaCoVAR(tt,i,j,qq,uu,6)=(outCoVAR(tt,i,j,qq,uu,6)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,6)=(outCoVAR(tt,j,i,qq,uu,6)-covarmedxj)/covarmedxj;
                    
                     % 17 Time-varying students't copula
                     rho17=dynamic_rho(tt,1,7,uu,i,j);
                     theta8tv2= dynamic_thetaALL(1,4,7,uu,i,j);
                    [fu,fval,exitflag]=fzero(@(u) copulacdf('t',[u QQ(qq)],rho17,theta8tv2)-u+alpha0-alpha0*beta0,[0 1]);
                    xi=tinv(fu, uandlamdai);
                    xj=tinv(fu, uandlamdaj);
                    outCoVAR(tt,i,j,qq,uu,7)=muhat(tt,i) + sqrt(hhat_opt(tt,i))*xi;  %%j对i的风险溢出
                    outCoVAR(tt,j,i,qq,uu,7)=muhat(tt,j) + sqrt(hhat_opt(tt,j))*xj;  %%i对j的风险溢出
                    
                    [f_medu,fval,exitflag]=fzero(@(u) copulacdf('t',[u 0.5],rho17,theta8tv2)-u+0.5-0.5*beta0,[0 1]);
                    medxi=tinv(f_medu, uandlamdai);
                    medxj=tinv(f_medu, uandlamdaj); 
                    covarmedxi=muhat(tt,i) + sqrt(hhat_opt(tt,i))*medxi;
                    covarmedxj=muhat(tt,j) + sqrt(hhat_opt(tt,j))*medxj;
                    outDeltaCoVAR(tt,i,j,qq,uu,7)=(outCoVAR(tt,i,j,qq,uu,7)-covarmedxi)/covarmedxi;
                    outDeltaCoVAR(tt,j,i,qq,uu,7)=(outCoVAR(tt,j,i,qq,uu,7)-covarmedxj)/covarmedxj;
                 end
            end
        end
        [i j toc]
    end
end


%%%%%%%%中国石油对石油期货的风险溢出
figure(5)
plot(datet,outVAR(:,1,2,2),'LineWidth',2)    %%5%分位数
hold on
plot(datet,outVAR(:,1,4,2),'LineWidth',2)    %%5%分位数
hold on

plot(datet,outCoVAR(:,1,2,2,2,5),'LineWidth',2)    %%5%分位数
hold on
plot(datet,outCoVAR(:,1,2,4,2,5),'LineWidth',2) 

legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'Best');grid on;



%%%%%%%%石油期货对中国石油的风险溢出
figure(6)
plot(datet,outVAR(:,2,2,2),'LineWidth',2)    %%5%分位数
hold on
plot(datet,outVAR(:,2,4,2),'LineWidth',2)    %%5%分位数
hold on

plot(datet,outCoVAR(:,2,1,2,2,5),'LineWidth',2)    %%5%分位数
hold on
plot(datet,outCoVAR(:,2,1,4,2,5),'LineWidth',2) 

legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'Best');grid on;

datetick('x','yy/mm/dd')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%表3 静态copula参数表

% 1. Normal Copula
thetaALL(:,1:5,1,2,1,2)
 % 2. Clayton's copula
 thetaALL(:,1:5,2,2,1,2)
 %3. Rotated Clayton copula
 thetaALL(:,1:5,3,2,1,2)
% 4. Gumbel copula
  thetaALL(:,1:5,4,2,1,2)  
   
% 5. Rotated Gumbel copula
  thetaALL(:,1:5,5,2,1,2)
% 6. Symmetrised Joe-Clayton copula
 thetaALL(:,1:5,6,2,1,2)

% 7. Student's t copula
 thetaALL(:,1:5,7,2,1,2)
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%表3 动态copula参数表

% 1. TVP-Normal Copula
dynamic_thetaALL(:,:,1,2,1,2)
 % 2. TVP-Clayton's copula
 dynamic_thetaALL(:,:,2,2,1,2)
 %3. TVP-Rotated Clayton copula
 dynamic_thetaALL(:,:,3,2,1,2)
% 4.TVP- Gumbel copula
 dynamic_thetaALL(:,:,4,2,1,2)  
   
% 5. TVP-Rotated Gumbel copula
  dynamic_thetaALL(:,:,5,2,1,2)
% 6. TVP-Symmetrised Joe-Clayton copula
 dynamic_thetaALL(:,:,6,2,1,2)

% 7. TVP-Student's t copula
 dynamic_thetaALL(:,:,7,2,1,2)
 
 
 
%%%%公司对行业的风险溢出
DescriptivetoC=nan(2,6,3);  %%1,2行为石油期货对中国石油的风险溢出统计
                           %%前3列为下5%分位数VaR,CoVaR,DeltaCoVaR，后三列为前3列为上5%分位数VaR,CoVaR,DeltaCoVaR
                           %%%3:表示3种分布 
for uu=2:3
    %%%%Downside，下5%处
    
    %%%石油期货对中国石油的风险溢出
    DescriptivetoC(1,1,uu)=mean(outVAR(:,2,2,uu));
    DescriptivetoC(2,1,uu)=std(outVAR(:,2,2,uu));
         
    numb1=dynamic_opt_copula(1,2,uu,1);  %%3为SPDB,24为银行业
    DescriptivetoC(1,2,uu)=mean(outCoVAR(:,2,1,2,uu,numb1));
    DescriptivetoC(2,2,uu)=std(outCoVAR(:,2,1,2,uu,numb1));
    DescriptivetoC(1,3,uu)=mean(outDeltaCoVAR(:,2,1,2,uu,numb1));
    DescriptivetoC(2,3,uu)=std(outDeltaCoVAR(:,2,1,2,uu,numb1));
         
      %%%%Upside,上95%
    DescriptivetoC(1,4,uu)=mean(outVAR(:,2,4,uu));
    DescriptivetoC(2,4,uu)=std(outVAR(:,2,4,uu));
         
    numb1=dynamic_opt_copula(1,2,uu,1);  %%3为SPDB,24为银行业
    DescriptivetoC(1,5,uu)=mean(outCoVAR(:,2,1,4,uu,numb1));
    DescriptivetoC(2,5,uu)=std(outCoVAR(:,2,1,4,uu,numb1));
    DescriptivetoC(1,6,uu)=mean(outDeltaCoVAR(:,2,1,4,uu,numb1));
    DescriptivetoC(2,6,uu)=std(outDeltaCoVAR(:,2,1,4,uu,numb1));
         
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%风险溢出检验，对应三个表，表1为风险溢出检验,表2为不同方向风险溢出大小检验


%%%%石油期货对中国石油的风险溢出检验

spillovertest1to=nan(2,4,3);   %%%第1列检验下5%分位数CoVaR<VaR,第2列检验上5%分位数CoVaR>VaR，
                                 %%%第3列检验CoVaR/VaR(Down)<CoVaR/VaR(Up),第4列检验DeltaCoVaR（Down）<DeltaCoVaR（Up）
                                 %%%奇数行为统计量，偶数行为统计量对应的P值        


for uu=2:3
    %%bank
     i=1;
     numb1=dynamic_opt_copula(1,2,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR(:,2,1,2,uu,numb1)', outVAR(:,2,2,uu)','Tail','larger');  %%%Downside
     spillovertest1to(2*i-1,1,uu)=k;
     spillovertest1to(2*i,1,uu)=p;
     
     [h p k]=kstest2(outCoVAR(:,2,1,4,uu,numb1)', outVAR(:,2,4,uu),'Tail','smaller');  %%%Upside
     spillovertest1to(2*i-1,2,uu)=k;
     spillovertest1to(2*i,2,uu)=p;
     
     [h p k]=kstest2(outCoVAR(:,2,1,2,uu,numb1)./outVAR(:,2,2,uu),outCoVAR(:,2,1,4,uu,numb1)./outVAR(:,2,4,uu),'Tail','smaller');
     spillovertest1to(2*i-1,3,uu)=k;
     spillovertest1to(2*i,3,uu)=p;
     
     [h p k]=kstest2(outDeltaCoVAR(:,2,1,2,uu,numb1),outDeltaCoVAR(:,2,1,4,uu,numb1),'Tail','smaller');
     spillovertest1to(2*i-1,4,uu)=k;
     spillovertest1to(2*i,4,uu)=p;
end

    
    











































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%以下代码暂时不用




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%CoVaR和DeltaCoVaR比较

%%%%%%%%%%%%%%%%%%
%%计算三个子行业中每个企业对行业的溢出效应，并根据ΔCoVaR识别溢出效应最大的企业

spill_DeltaCoVaR=nan(length(rets1),22,length(QQ),3);  %%其中前14个为银行业，15：19为多元金融业，20:22为保险业;3:对应对应分布类型
                                         %%spill_DeltaCoVaR为每个金融企业对各对应行业的溢出效应
for i=1:22
    for j=2:3   %%2表示偏t分布，3表示t分布
        for qq=1:length(QQ)
            if i<=14
                numb1=dynamic_opt_copula(i,24,j,1);  %%24为银行业
                spill_DeltaCoVaR(:,i,qq,j)=outDeltaCoVAR_ALL(:,24,i,qq,j,numb1);
            elseif i>=15 & i<=19
                numb1=dynamic_opt_copula(i,25,j,1);  %%25为多元金融业
                spill_DeltaCoVaR(:,i,qq,j)=outDeltaCoVAR_ALL(:,25,i,qq,j,numb1);
            elseif i>=20 & i<=22
                numb1=dynamic_opt_copula(i,26,j,1);  %%26为保险业
                spill_DeltaCoVaR(:,i,qq,j)=outDeltaCoVAR_ALL(:,26,i,qq,j,numb1);
            end
        end
    end
end





%%%对银行业溢出效应最大的
numb2=find(mean(spill_DeltaCoVaR(:,1:14,2,2))==max(mean(spill_DeltaCoVaR(:,1:14,2,2))))
sectors(numb2)

%%%对多元金融业溢出效应最大的
numb2=find(mean(spill_DeltaCoVaR(:,15:19,2,2))==max(mean(spill_DeltaCoVaR(:,15:19,2,2))))
sectors(14+numb2)


%%%对保险业溢出效应最大的
numb2=find(mean(spill_DeltaCoVaR(:,20:22,2,2))==max(mean(spill_DeltaCoVaR(:,20:22,2,2))))
sectors(19+numb2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%不同时间对行业溢出效应最大的公司

maxspill_company=nan(length(rets1),3,length(QQ),3);  %%第一个3中：第1列为溢出效应最大的银行，2列为多元金融，3列为保险；第2个3表示3中分布
maxspill_company2=cell(length(rets1),3,length(QQ),3);
for tt=1:length(rets1)
    for qq=[1:2 4:5]
        for uu=2:3
            numb1=find(spill_DeltaCoVaR(tt,1:14,qq,uu)==max(spill_DeltaCoVaR(tt,1:14,qq,uu)));
            maxspill_company(tt,1,qq,uu)=numb1;
            maxspill_company2(tt,1,qq,uu)=sectors(numb1);
            
            numb2=find(spill_DeltaCoVaR(tt,15:19,qq,uu)==max(spill_DeltaCoVaR(tt,15:19,qq,uu)));
            maxspill_company(tt,2,qq,uu)=14+numb2;
            maxspill_company2(tt,2,qq,uu)=sectors(14+numb2);
            
            numb3=find(spill_DeltaCoVaR(tt,20:22,qq,uu)==max(spill_DeltaCoVaR(tt,20:22,qq,uu)));
            maxspill_company(tt,3,qq,uu)=19+numb3(1);
            maxspill_company2(tt,3,qq,uu)=sectors(19+numb3(1));
        end
    end
end


%%%%%Figure 1.  Maximum Spillover Company Tracking Over Time.
opts.figPrefs = {'Units', 'Normalized', ...
    'Position', [0.125, 0.125, 0.75, 0.75], ...
    'Menubar', 'none'};
opts.titleSize = 14;
opts.labelSize = 12;

trackDates = datet(1:end);  %%%%%%%%修改，原命令为trackDates = dates(nWind+1:end);
figure(opts.figPrefs{:})
axLeft = subplot(1, 2, 1);

centralSectorIdx=maxspill_company(:,1,2,2);   %%%5%显著性水平下，偏t分布

plot(trackDates, centralSectorIdx, 'bo', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'y')
xlabel('Date', 'FontSize', opts.labelSize)
ylabel('Bank')
set(gca, 'YTick', 0:numel(sectors(1:14)), ...
    'YTickLabel', [{''}; sectors(1:14)'], ...  %%%%修改，原命令为[{''}; sectors]
    'FontSize', opts.labelSize, ...
    'YLim', [0, numel(sectors(1:14))+1], ...
    'XTickLabelRotation', 25)
grid

datetick('x','yy/mm/dd');
title('Maximum Spillover Bank Tracking Over Time')

%%%Maximum Spillover Diversified Financials
axTopRight = subplot(2, 2, 2);
centralSectorIdx=maxspill_company(:,2,2,2)-14;   %%%5%显著性水平下，偏t分布

plot(trackDates, centralSectorIdx, 'bo', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'y')
xlabel('Date', 'FontSize', opts.labelSize)
ylabel('Diversified Financials')
sectorsDF=sectors(15:19);    %%多元金融公司名称
set(gca, 'YTick', 0:numel(sectorsDF), ...
    'YTickLabel', [{''}; sectorsDF'], ...  %%%%修改，原命令为[{''}; sectors]
    'FontSize', opts.labelSize, ...
    'YLim', [0, numel(sectorsDF)+1], ...
    'XTickLabelRotation', 25)
grid

datetick('x','yyyy');
title('Maximum Spillover Diversified Financials Tracking Over Time')


%%%Maximum Spillover Insurance

axBottomRight = subplot(2, 2, 4);


centralSectorIdx=maxspill_company(:,3,2,2)-19;   %%%5%显著性水平下，偏t分布,多元金融

plot(trackDates, centralSectorIdx, 'bo', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'y')
xlabel('Date', 'FontSize', opts.labelSize)
ylabel('Insurance')
sectorsI=sectors(20:22);    %%多元金融公司名称
set(gca, 'YTick', 0:numel(sectorsI), ...
    'YTickLabel', [{''}; sectorsI'], ...  %%%%修改，原命令为[{''}; sectors]
    'FontSize', opts.labelSize, ...
    'YLim', [0, numel(sectorsI)+1], ...
    'XTickLabelRotation', 25)
grid

datetick('x','yyyy');
title('Maximum Spillover Insurance Tracking Over Time')








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%打开作图菜单栏,2018-5-11
shh = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','On')
set(gcf,'menubar','figure')
set(gcf,'CloseRequestFcn','closereq')
set(gcf,'DefaultLineClipping','Off')
set(0,'ShowHiddenHandles',shh)

%%title('')   %%%去掉标题

%%%%%%%%Figure 2. Distribution of Maximum Spillover Bank Tracking Over Time.
figure(opts.figPrefs{:})
axLeft = subplot(1, 2, 1);

centralSectorIdx=maxspill_company(:,1,2,2);   %%%5%显著性水平下，偏t分布
sectorsBank=sectors(1:14);    %%%%银行业中银行的名称


nodeCounts = accumarray(centralSectorIdx, ...
             ones(size(centralSectorIdx)), [numel(sectorsBank), 1], @sum, 0);
[nodeCounts, sortPos] = sort(nodeCounts);


barh(1:2:numel(nodeCounts), nodeCounts(1:2:end), 'FaceColor', 'b', ...
    'FaceAlpha', 0.5, 'BarWidth', 0.40, 'LineWidth', 1)
hold on
barh(2:2:numel(nodeCounts), nodeCounts(2:2:end), 'FaceColor', 'r', ...
    'FaceAlpha', 0.5, 'BarWidth', 0.40, 'LineWidth', 1)

set(gca, 'YTick', 1:numel(nodeCounts), ...
         'YTickLabel', sectorsBank(sortPos), ...
         'FontSize', opts.labelSize)
grid
title('Frequency of Maximum Spillover Bank', ...
    'FontSize', opts.titleSize)
nodeCountsPC = 100*nodeCounts/sum(nodeCounts);
textLabels = [repmat(' ', size(nodeCounts)), ...
    num2str(nodeCountsPC, '%.1f%%')];
text(nodeCounts, 1:numel(nodeCounts), textLabels, ...
    'Color', 'k', 'FontSize', opts.labelSize, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left')
%%%Maximum Spillover Diversified Financials
axTopRight = subplot(2, 2, 2);
centralSectorIdx=maxspill_company(:,2,2,2)-14;   %%%5%显著性水平下，偏t分布

sectorsDF=sectors(15:19);    %%%%银行业中银行的名称

nodeCounts = accumarray(centralSectorIdx, ...
             ones(size(centralSectorIdx)), [numel(sectorsDF), 1], @sum, 0);
[nodeCounts, sortPos] = sort(nodeCounts);


barh(1:2:numel(nodeCounts), nodeCounts(1:2:end), 'FaceColor', 'b', ...
    'FaceAlpha', 0.5, 'BarWidth', 0.40, 'LineWidth', 1)
hold on
barh(2:2:numel(nodeCounts), nodeCounts(2:2:end), 'FaceColor', 'r', ...
    'FaceAlpha', 0.5, 'BarWidth', 0.40, 'LineWidth', 1)

set(gca, 'YTick', 1:numel(nodeCounts), ...
         'YTickLabel', sectorsDF(sortPos), ...
         'FontSize', opts.labelSize)
grid
title('Frequency of Maximum Spillover Diversified Financials', ...
    'FontSize', opts.titleSize)
nodeCountsPC = 100*nodeCounts/sum(nodeCounts);
textLabels = [repmat(' ', size(nodeCounts)), ...
    num2str(nodeCountsPC, '%.1f%%')];
text(nodeCounts, 1:numel(nodeCounts), textLabels, ...
    'Color', 'k', 'FontSize', opts.labelSize, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left')

%%%Maximum Spillover Insurance

axBottomRight = subplot(2, 2, 4);
centralSectorIdx=maxspill_company(:,3,2,2)-19;   %%%5%显著性水平下，偏t分布,多元金融

sectorsI=sectors(20:22);    %%%%银行业中银行的名称

nodeCounts = accumarray(centralSectorIdx, ...
             ones(size(centralSectorIdx)), [numel(sectorsI), 1], @sum, 0);
[nodeCounts, sortPos] = sort(nodeCounts);


barh(1:2:numel(nodeCounts), nodeCounts(1:2:end), 'FaceColor', 'b', ...
    'FaceAlpha', 0.5, 'BarWidth', 0.40, 'LineWidth', 1)
hold on
barh(2:2:numel(nodeCounts), nodeCounts(2:2:end), 'FaceColor', 'r', ...
    'FaceAlpha', 0.5, 'BarWidth', 0.40, 'LineWidth', 1)

set(gca, 'YTick', 1:numel(nodeCounts), ...
         'YTickLabel', sectorsI(sortPos), ...
         'FontSize', opts.labelSize)
grid
title('Frequency of Maximum Spillover Insurance Tracking', ...
    'FontSize', opts.titleSize)
nodeCountsPC = 100*nodeCounts/sum(nodeCounts);
textLabels = [repmat(' ', size(nodeCounts)), ...
    num2str(nodeCountsPC, '%.1f%%')];
text(nodeCounts, 1:numel(nodeCounts), textLabels, ...
    'Color', 'k', 'FontSize', opts.labelSize, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%打开作图菜单栏,2018-5-11
shh = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','On')
set(gcf,'menubar','figure')
set(gcf,'CloseRequestFcn','closereq')
set(gcf,'DefaultLineClipping','Off')
set(0,'ShowHiddenHandles',shh)

%%title('')   %%%去掉标题


%%%%各公司处于某个时间点最大溢出效应统计
sectors

abc1=nan(14,2);
for i=1:14
    
    abc1(i,1)=i;
    abc1(i,2)=length(find(maxspill_company(:,1,2,2)==i));
end

abc2=nan(5,2);
for i=15:19
    
    abc2(i-14,1)=i;
    abc2(i-14,2)=length(find(maxspill_company(:,2,2,2)==i));
end


abc3=nan(3,2);
for i=20:22
    
    abc3(i-19,1)=i;
    abc3(i-19,2)=length(find(maxspill_company(:,3,2,2)==i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%风险溢出统计，只统计每个行业中对行业溢出效应最大的公司，其中银行业为SPDB,序号为3；多元金融业为PS，序号为16,；保险业为PAI，序号为20

%%%%公司对行业的风险溢出
DescriptivefromC=nan(3*2,6,3);  %%1,2行为SPDB对银行业的风险溢出统计，3,4行为PS对多元金融业；5,6行为PAI对保险业
                           %%前3列为下5%分位数VaR,CoVaR,DeltaCoVaR，后三列为前3列为上5%分位数VaR,CoVaR,DeltaCoVaR
                           %%%3:表示3种分布
for uu=2:3
         %%%%Downside
         
         %%bank
         DescriptivefromC(1,1,uu)=mean(outVAR(:,24,2,uu));
         DescriptivefromC(2,1,uu)=std(outVAR(:,24,2,uu));
         
         numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业
         DescriptivefromC(1,2,uu)=mean(outCoVAR_ALL(:,24,3,2,uu,numb1));
         DescriptivefromC(2,2,uu)=std(outCoVAR_ALL(:,24,3,2,uu,numb1));
         DescriptivefromC(1,3,uu)=mean(outDeltaCoVAR_ALL(:,24,3,2,uu,numb1));
         DescriptivefromC(2,3,uu)=std(outDeltaCoVAR_ALL(:,24,3,2,uu,numb1));
         
         %%Diversified Financials
         DescriptivefromC(3,1,uu)=mean(outVAR(:,25,2,uu));
         DescriptivefromC(4,1,uu)=std(outVAR(:,25,2,uu));
         
         numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业
         DescriptivefromC(3,2,uu)=mean(outCoVAR_ALL(:,25,16,2,uu,numb1));
         DescriptivefromC(4,2,uu)=std(outCoVAR_ALL(:,25,16,2,uu,numb1));
         DescriptivefromC(3,3,uu)=mean(outDeltaCoVAR_ALL(:,25,16,2,uu,numb1));
         DescriptivefromC(4,3,uu)=std(outDeltaCoVAR_ALL(:,25,16,2,uu,numb1));
         
         %%Insurance
         DescriptivefromC(5,1,uu)=mean(outVAR(:,26,2,uu));
         DescriptivefromC(6,1,uu)=std(outVAR(:,26,2,uu));
         
         numb1=dynamic_opt_copula(20,26,uu,1);  %%3为SPDB,24为银行业
         DescriptivefromC(5,2,uu)=mean(outCoVAR_ALL(:,26,20,2,uu,numb1));
         DescriptivefromC(6,2,uu)=std(outCoVAR_ALL(:,26,20,2,uu,numb1));
         DescriptivefromC(5,3,uu)=mean(outDeltaCoVAR_ALL(:,26,20,2,uu,numb1));
         DescriptivefromC(6,3,uu)=std(outDeltaCoVAR_ALL(:,26,20,2,uu,numb1));
         
         %%%%Upside
         
         %%bank
         DescriptivefromC(1,4,uu)=mean(outVAR(:,24,4,uu));
         DescriptivefromC(2,4,uu)=std(outVAR(:,24,4,uu));
         
         numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业
         DescriptivefromC(1,5,uu)=mean(outCoVAR_ALL(:,24,3,4,uu,numb1));
         DescriptivefromC(2,5,uu)=std(outCoVAR_ALL(:,24,3,4,uu,numb1));
         DescriptivefromC(1,6,uu)=mean(outDeltaCoVAR_ALL(:,24,3,4,uu,numb1));
         DescriptivefromC(2,6,uu)=std(outDeltaCoVAR_ALL(:,24,3,4,uu,numb1));
         
         %%Diversified Financials
         DescriptivefromC(3,4,uu)=mean(outVAR(:,25,4,uu));
         DescriptivefromC(4,4,uu)=std(outVAR(:,25,4,uu));
         
         numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业
         DescriptivefromC(3,5,uu)=mean(outCoVAR_ALL(:,25,16,4,uu,numb1));
         DescriptivefromC(4,5,uu)=std(outCoVAR_ALL(:,25,16,4,uu,numb1));
         DescriptivefromC(3,6,uu)=mean(outDeltaCoVAR_ALL(:,25,16,4,uu,numb1));
         DescriptivefromC(4,6,uu)=std(outDeltaCoVAR_ALL(:,25,16,4,uu,numb1));
         
         %%Insurance
         DescriptivefromC(5,4,uu)=mean(outVAR(:,26,4,uu));
         DescriptivefromC(6,4,uu)=std(outVAR(:,26,4,uu));
         
         numb1=dynamic_opt_copula(20,26,uu,1);  %%3为SPDB,24为银行业
         DescriptivefromC(5,5,uu)=mean(outCoVAR_ALL(:,26,20,4,uu,numb1));
         DescriptivefromC(6,5,uu)=std(outCoVAR_ALL(:,26,20,4,uu,numb1));
         DescriptivefromC(5,6,uu)=mean(outDeltaCoVAR_ALL(:,26,20,4,uu,numb1));
         DescriptivefromC(6,6,uu)=std(outDeltaCoVAR_ALL(:,26,20,4,uu,numb1));
end


%%%%行业对公司的风险溢出
DescriptivetoC=nan(3*2,6,3);  %%1,2行为SPDB对银行业的风险溢出统计，3,4行为PS对多元金融业；5,6行为PAI对保险业
                           %%前3列为下5%分位数VaR,CoVaR,DeltaCoVaR，后三列为前3列为上5%分位数VaR,CoVaR,DeltaCoVaR
                           %%%3:表示3种分布
for uu=2:3
         %%%%Downside
         
         %%bank
         DescriptivetoC(1,1,uu)=mean(outVAR(:,3,2,uu));
         DescriptivetoC(2,1,uu)=std(outVAR(:,3,2,uu));
         
         numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业
         DescriptivetoC(1,2,uu)=mean(outCoVAR_ALL(:,3,24,2,uu,numb1));
         DescriptivetoC(2,2,uu)=std(outCoVAR_ALL(:,3,24,2,uu,numb1));
         DescriptivetoC(1,3,uu)=mean(outDeltaCoVAR_ALL(:,3,24,2,uu,numb1));
         DescriptivetoC(2,3,uu)=std(outDeltaCoVAR_ALL(:,3,24,2,uu,numb1));
         
         %%Diversified Financials
         DescriptivetoC(3,1,uu)=mean(outVAR(:,16,2,uu));
         DescriptivetoC(4,1,uu)=std(outVAR(:,16,2,uu));
         
         numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业
         DescriptivetoC(3,2,uu)=mean(outCoVAR_ALL(:,16,25,2,uu,numb1));
         DescriptivetoC(4,2,uu)=std(outCoVAR_ALL(:,16,25,2,uu,numb1));
         DescriptivetoC(3,3,uu)=mean(outDeltaCoVAR_ALL(:,16,25,2,uu,numb1));
         DescriptivetoC(4,3,uu)=std(outDeltaCoVAR_ALL(:,16,25,2,uu,numb1));
         
         %%Insurance
         DescriptivetoC(5,1,uu)=mean(outVAR(:,20,2,uu));
         DescriptivetoC(6,1,uu)=std(outVAR(:,20,2,uu));
         
         numb1=dynamic_opt_copula(20,26,uu,1);  %%3为SPDB,24为银行业
         DescriptivetoC(5,2,uu)=mean(outCoVAR_ALL(:,20,26,2,uu,numb1));
         DescriptivetoC(6,2,uu)=std(outCoVAR_ALL(:,20,26,2,uu,numb1));
         DescriptivetoC(5,3,uu)=mean(outDeltaCoVAR_ALL(:,20,26,2,uu,numb1));
         DescriptivetoC(6,3,uu)=std(outDeltaCoVAR_ALL(:,20,26,2,uu,numb1));
         
         %%%%Upside
         
         %%bank
         DescriptivetoC(1,4,uu)=mean(outVAR(:,3,4,uu));
         DescriptivetoC(2,4,uu)=std(outVAR(:,3,4,uu));
         
         numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业
         DescriptivetoC(1,5,uu)=mean(outCoVAR_ALL(:,3,24,4,uu,numb1));
         DescriptivetoC(2,5,uu)=std(outCoVAR_ALL(:,3,24,4,uu,numb1));
         DescriptivetoC(1,6,uu)=mean(outDeltaCoVAR_ALL(:,3,24,4,uu,numb1));
         DescriptivetoC(2,6,uu)=std(outDeltaCoVAR_ALL(:,3,24,4,uu,numb1));
         
         %%Diversified Financials
         DescriptivetoC(3,4,uu)=mean(outVAR(:,16,4,uu));
         DescriptivetoC(4,4,uu)=std(outVAR(:,16,4,uu));
         
         numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业
         DescriptivetoC(3,5,uu)=mean(outCoVAR_ALL(:,16,25,4,uu,numb1));
         DescriptivetoC(4,5,uu)=std(outCoVAR_ALL(:,16,25,4,uu,numb1));
         DescriptivetoC(3,6,uu)=mean(outDeltaCoVAR_ALL(:,16,25,4,uu,numb1));
         DescriptivetoC(4,6,uu)=std(outDeltaCoVAR_ALL(:,16,25,4,uu,numb1));
         
         %%Insurance
         DescriptivetoC(5,4,uu)=mean(outVAR(:,20,4,uu));
         DescriptivetoC(6,4,uu)=std(outVAR(:,20,4,uu));
         
         numb1=dynamic_opt_copula(20,26,uu,1);  %%3为SPDB,24为银行业
         DescriptivetoC(5,5,uu)=mean(outCoVAR_ALL(:,20,26,4,uu,numb1));
         DescriptivetoC(6,5,uu)=std(outCoVAR_ALL(:,20,26,4,uu,numb1));
         DescriptivetoC(5,6,uu)=mean(outDeltaCoVAR_ALL(:,20,26,4,uu,numb1));
         DescriptivetoC(6,6,uu)=std(outDeltaCoVAR_ALL(:,20,26,4,uu,numb1));
end


Descriptivefromandto=[DescriptivefromC(:,:,2)
    DescriptivetoC(:,:,2)];

xlswrite('风险溢出描述性统计1.xlsx',Descriptivefromandto)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%风险溢出检验，对应三个表，表1为风险溢出检验,表2为不同方向风险溢出大小检验


%%%%公司对行业的风险溢出
spillovertest1from=nan(3*2,4,3);   %%%第1列检验下5%分位数CoVaR<VaR,第2列检验上5%分位数CoVaR>VaR，
                                 %%%第3列检验CoVaR/VaR(Down)<CoVaR/VaR(Up),第4列检验DeltaCoVaR（Down）<DeltaCoVaR（Up）
                                 %%%奇数行为统计量，偶数行为统计量对应的P值  
                                 %%%3表示分布

for uu=2:3
     %%bank
     i=1;
     numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR_ALL(:,24,3,2,uu,numb1)', outVAR(:,24,2,uu)','Tail','larger');  %%%Downside
     spillovertest1from(2*i-1,1,uu)=k;
     spillovertest1from(2*i,1,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,24,3,4,uu,numb1)', outVAR(:,24,4,uu),'Tail','smaller');  %%%Upside
     spillovertest1from(2*i-1,2,uu)=k;
     spillovertest1from(2*i,2,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,24,3,2,uu,numb1)./outVAR(:,24,2,uu),outCoVAR_ALL(:,24,3,4,uu,numb1)./outVAR(:,24,4,uu),'Tail','larger');
     spillovertest1from(2*i-1,3,uu)=k;
     spillovertest1from(2*i,3,uu)=p;
     
     [h p k]=kstest2(outDeltaCoVAR_ALL(:,24,3,2,uu,numb1),outDeltaCoVAR_ALL(:,24,3,4,uu,numb1),'Tail','larger');
     spillovertest1from(2*i-1,4,uu)=k;
     spillovertest1from(2*i,4,uu)=p;
     
     %%Diversified Financials
     i=2;
     numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR_ALL(:,25,16,2,uu,numb1), outVAR(:,25,2,uu),'Tail','larger');  %%%Downside
     spillovertest1from(2*i-1,1,uu)=k;
     spillovertest1from(2*i,1,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,25,16,4,uu,numb1)', outVAR(:,25,4,uu),'Tail','smaller');  %%%Upside
     spillovertest1from(2*i-1,2,uu)=k;
     spillovertest1from(2*i,2,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,25,16,2,uu,numb1)./outVAR(:,25,2,uu),outCoVAR_ALL(:,25,16,4,uu,numb1)./outVAR(:,25,4,uu),'Tail','larger');
     spillovertest1from(2*i-1,3,uu)=k;
     spillovertest1from(2*i,3,uu)=p;
     
     [h p k]=kstest2(outDeltaCoVAR_ALL(:,25,16,2,uu,numb1),outDeltaCoVAR_ALL(:,25,16,4,uu,numb1),'Tail','larger');
     spillovertest1from(2*i-1,4,uu)=k;
     spillovertest1from(2*i,4,uu)=p;
     
     %%Insurance
      i=3;
     numb1=dynamic_opt_copula(20,26,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR_ALL(:,26,20,2,uu,numb1), outVAR(:,26,2,uu),'Tail','larger');  %%%Downside
     spillovertest1from(2*i-1,1,uu)=k;
     spillovertest1from(2*i,1,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,26,20,4,uu,numb1)', outVAR(:,26,4,uu),'Tail','smaller');  %%%Upside
     spillovertest1from(2*i-1,2,uu)=k;
     spillovertest1from(2*i,2,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,26,20,2,uu,numb1)./outVAR(:,26,2,uu),outCoVAR_ALL(:,26,20,4,uu,numb1)./outVAR(:,26,4,uu),'Tail','larger');
     spillovertest1from(2*i-1,3,uu)=k;
     spillovertest1from(2*i,3,uu)=p;
     
     [h p k]=kstest2(outDeltaCoVAR_ALL(:,26,20,2,uu,numb1),outDeltaCoVAR_ALL(:,26,20,4,uu,numb1),'Tail','larger');
     spillovertest1from(2*i-1,4,uu)=k;
     spillovertest1from(2*i,4,uu)=p;
     
end

%%%%行业对公司的风险溢出
spillovertest1to=nan(3*2,4,3);   %%%第1列检验下5%分位数CoVaR<VaR,第2列检验上5%分位数CoVaR>VaR，
                                 %%%第3列检验CoVaR/VaR(Down)<CoVaR/VaR(Up),第4列检验DeltaCoVaR（Down）<DeltaCoVaR（Up）
                                 %%%奇数行为统计量，偶数行为统计量对应的P值        

for uu=2:3
     %%bank
     i=1;
     numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR_ALL(:,3,24,2,uu,numb1)', outVAR(:,3,2,uu)','Tail','larger');  %%%Downside
     spillovertest1to(2*i-1,1,uu)=k;
     spillovertest1to(2*i,1,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,3,24,4,uu,numb1)', outVAR(:,3,4,uu),'Tail','smaller');  %%%Upside
     spillovertest1to(2*i-1,2,uu)=k;
     spillovertest1to(2*i,2,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,3,24,2,uu,numb1)./outVAR(:,3,2,uu),outCoVAR_ALL(:,3,24,4,uu,numb1)./outVAR(:,3,4,uu),'Tail','larger');
     spillovertest1to(2*i-1,3,uu)=k;
     spillovertest1to(2*i,3,uu)=p;
     
     [h p k]=kstest2(outDeltaCoVAR_ALL(:,3,24,2,uu,numb1),outDeltaCoVAR_ALL(:,3,24,4,uu,numb1),'Tail','larger');
     spillovertest1to(2*i-1,4,uu)=k;
     spillovertest1to(2*i,4,uu)=p;
     
     %%Diversified Financials
     i=2;
     numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR_ALL(:,16,25,2,uu,numb1), outVAR(:,16,2,uu),'Tail','larger');  %%%Downside
     spillovertest1to(2*i-1,1,uu)=k;
     spillovertest1to(2*i,1,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,16,25,4,uu,numb1)', outVAR(:,16,4,uu),'Tail','smaller');  %%%Upside
     spillovertest1to(2*i-1,2,uu)=k;
     spillovertest1to(2*i,2,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,16,25,2,uu,numb1)./outVAR(:,16,2,uu),outCoVAR_ALL(:,16,25,4,uu,numb1)./outVAR(:,16,4,uu),'Tail','larger');
     spillovertest1to(2*i-1,3,uu)=k;
     spillovertest1to(2*i,3,uu)=p;
     
     [h p k]=kstest2(outDeltaCoVAR_ALL(:,16,25,2,uu,numb1),outDeltaCoVAR_ALL(:,16,25,4,uu,numb1),'Tail','larger');
     spillovertest1to(2*i-1,4,uu)=k;
     spillovertest1to(2*i,4,uu)=p;
     
     %%Insurance
      i=3;
     numb1=dynamic_opt_copula(20,26,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR_ALL(:,20,26,2,uu,numb1), outVAR(:,20,2,uu),'Tail','larger');  %%%Downside
     spillovertest1to(2*i-1,1,uu)=k;
     spillovertest1to(2*i,1,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,20,26,4,uu,numb1)', outVAR(:,20,4,uu),'Tail','smaller');  %%%Upside
     spillovertest1to(2*i-1,2,uu)=k;
     spillovertest1to(2*i,2,uu)=p;
     
     [h p k]=kstest2(outCoVAR_ALL(:,20,26,2,uu,numb1)./outVAR(:,20,2,uu),outCoVAR_ALL(:,20,26,4,uu,numb1)./outVAR(:,20,4,uu),'Tail','larger');
     spillovertest1to(2*i-1,3,uu)=k;
     spillovertest1to(2*i,3,uu)=p;
     
     [h p k]=kstest2(outDeltaCoVAR_ALL(:,20,26,2,uu,numb1),outDeltaCoVAR_ALL(:,20,26,4,uu,numb1),'Tail','larger');
     spillovertest1to(2*i-1,4,uu)=k;
     spillovertest1to(2*i,4,uu)=p;
     
end


spillovertest1fromandto=[spillovertest1from(:,:,2)
    spillovertest1to(:,:,2)];
xlswrite('风险溢出检验1.xlsx',spillovertest1fromandto) 


%%%%不同溢出方向比较
spillovertest2=nan(2*3,2,3);    %%%%判断from IBC和 to IBC的差别
                               %%第一列为Downside，第二列为Upside,3:3种分布

 for uu=2:3
     numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业
     [h p k]=kstest2(outCoVAR_ALL(:,24,3,2,uu,numb1), outCoVAR_ALL(:,3,24,2,uu,numb1),'Tail','smaller');  %%%Downside
     spillovertest2(2*i-1,1)=k;
     spillovertest2(2*i,1)=p;
     [h p k]=kstest2(outCoVAR_ALL(:,24,3,4,uu,numb1), outCoVAR_ALL(:,3,24,4,uu,numb1),'Tail','smaller');  %%%Downside
     spillovertest2(2*i-1,1)=k;
     spillovertest2(2*i,1)=p;
                               
     
     [h p k]=kstest2(outCoVAR_ALL(:,24,3,2,uu,numb1),outCoVAR_ALL(:,24,3,2,uu,numb1),'Tail','smaller');
     spillovertest3(2*i-1,1)=k;
     spillovertest3(2*i,1)=p;
     [h p k]=kstest2(outDeltaCoVAR_tfrom8(i,:,4,2,3)',outDeltaCoVAR_tto8(i,:,4,2,3)','Tail','smaller');
     spillovertest3(2*i-1,2)=k;
     spillovertest3(2*i,2)=p;
 end
 
 
 
xlswrite('同行业不同方向风险溢出比较.xlsx',spillovertest2)   
                              


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%风险溢出动态图、VaR和CoVaR


%%%%公司对行业的风险溢出

%%bank
uu=2;
numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业

figure(1),plot(datet,outVAR(:,24,2,uu),'b-','LineWidth',2);hold on;
plot(datet,outVAR(:,24,4,uu),'r--','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,24,3,2,uu,numb1),'k:','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,24,3,4,uu,numb1),'g-.','LineWidth',2);hold on;
legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'NorthEast');grid on;
title('Banking', 'FontSize', 12)
datetick('x','yyyy');


%%Diversified Financials
uu=2;
numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业


figure(2),plot(datet,outVAR(:,25,2,uu),'b-','LineWidth',2);hold on;
plot(datet,outVAR(:,25,4,uu),'r--','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,25,16,2,uu,numb1),'k:','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,25,16,4,uu,numb1),'g-.','LineWidth',2);hold on;
legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'NorthEast');grid on;
title('Diversified Financials', 'FontSize', 12)
datetick('x','yyyy');

 %%Insurance
uu=2;
numb1=dynamic_opt_copula(20,26,uu,1);  %%


figure(3),plot(datet,outVAR(:,26,2,uu),'b-','LineWidth',2);hold on;
plot(datet,outVAR(:,26,4,uu),'r--','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,26,20,2,uu,numb1),'k:','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,26,20,4,uu,numb1),'g-.','LineWidth',2);hold on;
legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'NorthEast');grid on;
title('Insurance', 'FontSize', 12)
datetick('x','yyyy');


%%%%行业对公司的风险溢出

%%bank
uu=2;
numb1=dynamic_opt_copula(3,24,uu,1);  %%3为SPDB,24为银行业

figure(1),plot(datet,outVAR(:,3,2,uu),'b-','LineWidth',2);hold on;
plot(datet,outVAR(:,3,4,uu),'r--','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,3,24,2,uu,numb1),'k:','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,3,24,4,uu,numb1),'g-.','LineWidth',2);hold on;
legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'NorthEast');grid on;
title('Banking', 'FontSize', 12)
datetick('x','yyyy');


%%Diversified Financials
uu=2;
numb1=dynamic_opt_copula(16,25,uu,1);  %%3为SPDB,24为银行业


figure(2),plot(datet,outVAR(:,16,2,uu),'b-','LineWidth',2);hold on;
plot(datet,outVAR(:,16,4,uu),'r--','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,16,25,2,uu,numb1),'k:','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,16,25,4,uu,numb1),'g-.','LineWidth',2);hold on;
legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'NorthEast');grid on;
title('Diversified Financials', 'FontSize', 12)
datetick('x','yyyy');

 %%Insurance
uu=2;
numb1=dynamic_opt_copula(20,26,uu,1);  %%


figure(3),plot(datet,outVAR(:,20,2,uu),'b-','LineWidth',2);hold on;
plot(datet,outVAR(:,20,4,uu),'r--','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,20,26,2,uu,numb1),'k:','LineWidth',2);hold on;
plot(datet,outCoVAR_ALL(:,20,26,4,uu,numb1),'g-.','LineWidth',2);hold on;
legend('VaR(D)','VaR(U)','CoVaR(D)','CoVaR(U)','Location', 'NorthEast');grid on;
title('Insurance', 'FontSize', 12)
datetick('x','yyyy');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%风险溢出动态图、DeltaCoVaR


%%%%公司对行业的风险溢出

%%bank
uu=2;
numb1=dynamic_opt_copula(3,24,uu,1);  %%

figure(4),plot(datet(22:end),outDeltaCoVAR_ALL(22:end,24,3,2,uu,numb1),'b-','LineWidth',1);hold on;  %%5%分位数
plot(datet(25:end),outDeltaCoVAR_ALL(25:end,24,3,4,uu,numb1),'r--','LineWidth',1);  %%95%分位数

legend('DeltaCoVaR(D)','DeltaCoVaR(U)','Location', 'NorthEast');grid on;
title('Banking', 'FontSize', 12)
datetick('x','yyyy');

%%Diversified Financials
uu=2;
numb1=dynamic_opt_copula(16,25,uu,1);  %%


figure(5),plot(datet(16:end),outDeltaCoVAR_ALL(16:end,25,16,2,uu,numb1),'b-','LineWidth',1);hold on;  %%5%分位数
plot(datet(16:end),outDeltaCoVAR_ALL(16:end,25,16,4,uu,numb1),'r--','LineWidth',1);  %%95%分位数

legend('DeltaCoVaR(D)','DeltaCoVaR(U)','Location', 'NorthEast');grid on;
title('Diversified Financials', 'FontSize', 12)
datetick('x','yyyy');

%%Insurance
uu=2;
numb1=dynamic_opt_copula(20,26,uu,1);  %%


figure(6),plot(datet(13:end),outDeltaCoVAR_ALL(13:end,26,20,2,uu,numb1),'b-','LineWidth',1);hold on;  %%5%分位数
plot(datet(13:end),outDeltaCoVAR_ALL(13:end,26,20,4,uu,numb1),'r--','LineWidth',1);  %%95%分位数

legend('DeltaCoVaR(D)','DeltaCoVaR(U)','Location', 'NorthEast');grid on;
title('Insurance', 'FontSize', 12)
datetick('x','yyyy');
                                 
                                 
                                 
                                                           
%%%%行业对公司的风险溢出
 %%bank
uu=2;
numb1=dynamic_opt_copula(3,24,uu,1);  %%

figure(7),plot(datet(11:end),outDeltaCoVAR_ALL(11:end,3,24,2,uu,numb1),'b-','LineWidth',1);hold on;  %%5%分位数
plot(datet(11:end),outDeltaCoVAR_ALL(11:end,3,24,4,uu,numb1),'r--','LineWidth',1);  %%95%分位数

legend('DeltaCoVaR(D)','DeltaCoVaR(U)','Location', 'NorthEast');grid on;
title('Banking', 'FontSize', 12)
datetick('x','yyyy');

%%Diversified Financials
uu=2;
numb1=dynamic_opt_copula(16,25,uu,1);  %%


figure(8),plot(datet(11:end),outDeltaCoVAR_ALL(11:end,16,25,2,uu,numb1),'b-','LineWidth',1);hold on;  %%5%分位数
plot(datet(11:end),outDeltaCoVAR_ALL(11:end,16,25,4,uu,numb1),'r--','LineWidth',1);  %%95%分位数

legend('DeltaCoVaR(D)','DeltaCoVaR(U)','Location', 'NorthEast');grid on;
title('Diversified Financials', 'FontSize', 12)
datetick('x','yyyy');

%%Insurance
uu=2;
numb1=dynamic_opt_copula(20,26,uu,1);  %%


figure(9),plot(datet(15:end),outDeltaCoVAR_ALL(15:end,20,26,2,uu,numb1),'b-','LineWidth',1);hold on;  %%5%分位数
plot(datet(11:end),outDeltaCoVAR_ALL(11:end,20,26,4,uu,numb1),'r--','LineWidth',1);  %%95%分位数

legend('DeltaCoVaR(D)','DeltaCoVaR(U)','Location', 'NorthEast');grid on;
title('Insurance', 'FontSize', 12)
datetick('x','yyyy');
                             

