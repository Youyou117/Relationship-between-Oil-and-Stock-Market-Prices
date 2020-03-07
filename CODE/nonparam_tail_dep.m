function [tauLU,tauLUQ, tauBOOT, QQb] = nonparam_tail_dep(data,QQ,bootreps,alpha,fig);
% function [tauL, tauU, tauLQ, tauUQ, tauBOOT, QQb] = nonparam_tail_dep(data,QQ,bootreps,alpha,fig);
%
%  Function to compute the nonparametric estimators of tail dependence considered in Frahm, Junker and Schmidt (2005, IME)
%
%  INPUTS:  data, a Tx2 matrix of data
%           QQ, a Qx1 vector of cut-offs to use (default = (1/T:1/T:0.1)')
%           bootreps, a scalar, the number of bootstrap replications to use to get a bootstrap confidence interval (default=0)
%           alpha, a scalar, the signif level for the bootstrap confidence interval (default = 0.05)
%           fig, a scalar, =1 if want a plot of estimated tail dep as a function of cut-off q, =0 else (default=0)
%
%  OUTPUTS: tauLU, a 3x2 vector, the estimated lower and upper tail dep coefficients, using the "log" estimator and the optimal threshold selection rule in Frahm et al.
%           tauLUQ, a Qx2 vector, the estimated lower and upper tail dep coefficients for each of the cut-offs
%           tauBOOT, a bootrepsx2 matrix, the bootstrap estimates of upper and lower tail dep
%           QQb, a Qx1 vector, the vector of cut-offs used 
%
%  Andrew Patton
%
%  12 April 2012

% Note: I implement the "log" estimator from Frahm, et al. below. In simulations and empirical work I got similar estimates from the "SEC" estimator
% in that paper. The "CFG" estimator in that paper requires taking "block maxima" and can give quite different estimates without that step.

%%
T = size(data,1);
b = floor(0.005*T);  % smoothing parameter (p94)
b = min(b,10);  % AJP: making sure we don't smooth over too many estimates.

data = empiricalCDF(data);  % transforming data to have Unif(0,1) margins

if nargin<2 || isempty(QQ);
    QQ = linspace(2/T,0.1,100+b)';
end
Q = length(QQ);
b = min(b,floor(Q/3));

if nargin<3 || isempty(bootreps)
    bootreps=0;
end
if nargin<4 || isempty(alpha)
    alpha=0.05;
end
if nargin<5 || isempty(fig)
    fig=0;
end

[tauLU,tauLUQ] = nonparam_tail_dep_calc(data,QQ,b);  % point estimates of lower and upper tail dep

%%
if bootreps>0
    
    bootdates = randi(T,T,bootreps);  % a set of dates to use in bootstrap (data are assumed iid so just shuffle dates without using blocks
    
    tauBOOT = nan(bootreps,2);  % [bootreps] ; [log, SEC estimators] ; [lower, upper tail]
    for bb=1:bootreps
        tauLUa = nonparam_tail_dep_calc(data(bootdates(:,bb),:),QQ,b);
        tauBOOT(bb,:) = tauLUa;
    end
    
    for tt=1:2;
        tauLU(2,tt) = quantile(tauBOOT(:,tt),alpha/2);
        tauLU(3,tt) = quantile(tauBOOT(:,tt),1-alpha/2);
    end
end

%%
if fig==1 
    figure(101);
    subplot(2,1,1);plot(QQ,tauLUQ(:,1),'bo-');hold on;
    plot(QQ,tauLU(1,1)*ones(Q,1),'b--');
    if bootreps>0;
        jbfill(QQ',tauLU(2,1)*ones(1,Q),tauLU(3,1)*ones(1,Q),'b','b',[],0.2);hold on;
    end
    axis([0,0.1,0,1]),legend('Estimate for various q','Point estimate','90% conf. int.',1),title('Lower tail dependence'),xlabel('Cut-off threshold (q)'),ylabel('Tail dep');hold off;
    
    subplot(2,1,2);plot(QQ,tauLUQ(:,1),'rp-');hold on;
    plot(QQ,tauLU(1,2)*ones(Q,1),'r--');
    if bootreps>0;
        jbfill(QQ',tauLU(2,2)*ones(1,Q),tauLU(3,2)*ones(1,Q),'r','r',[],0.2);hold on;
    end
    axis([0,0.1,0,1]),legend('Estimate for various q','Point estimate','90% conf. int.',1),title('Upper tail dependence'),xlabel('Cut-off threshold (q)'),ylabel('Tail dep');hold off;
end
        

%%
% helper function: computes the quantile dependence
    function [tauLUb,tauLUQb] = nonparam_tail_dep_calc(dataU,QQ,b)
        
        Q = length(QQ);
        tauLUb = nan(1,2);  % rows: [estimate; conf int LB ; conf int UB ]   cols: [lower tail, upper tail]
        tauLUQb = nan(Q,2);
        
        dataUrot = 1-dataU;
        % estimating the tail dep coefficients for a range of values of Q
        temp = empirical_copula(dataUrot(:,1:2),[1-QQ,1-QQ]);  % empirical copula evaluated at (q,q)
        tauLUQb(:,1) = 2-log(temp)./log(1-QQ);   % "log" estimator in Frahm et al. (2005, IME) p84
        tauLUQb(:,1) = min(max(tauLUQb(:,1),0),1);  % making sure estimates are inside [0,1]
       
        temp = empirical_copula(dataU(:,1:2),[1-QQ,1-QQ]);  % empirical copula evaluated at (q,q)
        tauLUQb(:,2) = 2-log(temp)./log(1-QQ);   % "log" estimator in Frahm et al. (2005, IME) p84
        tauLUQb(:,2) = min(max(tauLUQb(:,2),0),1);  % making sure estimates are inside [0,1]

        % now smoothing these across b choices of cut-off q, as in Frahm et al.:
        for tt=1:2;
            temp = [tauLUQb(:,tt),mlag(tauLUQb(:,tt),b-1)];
            tauLUQa(:,tt) = mean(temp(b+1:end,:)')';
        end
        
        m = floor(sqrt(T-2*b));  % desired length of plateau (p95)
        m = min(m,20);  % AJP: don't want the plateau to be too large
        m = min(m,floor(Q/3));  % can't smooth acros avalues of Q unless we have quite a few.
        for tt=1:2;  % now looking for "plateau" in estimated tail dep (ie, a flat spot as a function of q)
            SIG = std(tauLUQa(:,tt));
            temp = [tauLUQa(:,tt),mlag(tauLUQa(:,tt),m-1)];
            temp = temp(m+1:end,:);
            tempABS = sum(abs(diff(temp')')')';  % first eqn of p95
            kstar = find(tempABS<=2*SIG);
            if length(kstar)>1
                kstar = min(kstar);
            elseif length(kstar)==0  % ie, this condition is never satisfied)
                kstar = 1;  % then I will choose the first plateau. (Frahm et al. suggest setting tail dep to zero in this case, but I think that is too harsh.)
            end
            tauLUb(tt) = tauLUQa(kstar,tt);
        end
        
        
    end

end