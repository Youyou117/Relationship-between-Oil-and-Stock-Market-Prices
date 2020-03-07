function [out1,qc1,teststat,pval] = exceedence_correl(X,Y,qc,indic)
% Computes the 'exceedence correlations' discussed in
% Ang and Chen, Journal of Financial Economics, 2002.
%
% Corr[X,Y|X>v,Y>v] if v>0 and X<v,Y<v if v<0
%
% Testing the significance of the differences using the test of Hong, Tu and Zhou (2003).
% 
% INPUTS:	X, a Tx1 vector of data
%				Y, a Tx1 vector of data
%				qc, a (k/2)x1 vector of quantiles or standard deviations to estimate the correl. must be increasing, and only the
%		    		"upper" part (ie, >=0.5 for q and >=0 for c)
%               indic, a scalar, equals 1 if using quantiles to determine cut-offs, =0 if using standard deviations
%
% OUTPUT:	out1, a (k+1)x1 vector of the measure at each quantile
%
% Saturday, 13 December, 2003
%
% Andrew Patton
%
% see also: quantiledep.m

if nargin<4
    indic=1;    % using quantiles to determine cut-offs
end
if nargin<3
    qc=0.5;     % just looking at above or below median returns
end
T = size(X,1);
k = size(qc,1);
if indic==1
    out1 = nines(2*k,1);		
    
    % firstly: ripping through to make sure I have enough data in each tail
    qc2 = qc;
    for jj=1:length(qc);
        temp1 = find( (X<=quantile(X,1-qc(jj))).*(Y<=quantile(Y,1-qc(jj))) );
        temp2 = find( (X>=quantile(X,qc(jj))).*(Y>=quantile(Y,qc(jj))) );
        if ((length(temp1)<3)+(length(temp2)<3))>0  % then not enough data to compute a correlation in at least one tail 
            qc2 = setdiff(qc2,qc(jj));
        end
    end
    
    % now I am sure that each cut-off has enough data
    qc2 = sortrows(unique([qc2;1-qc2]));   % sorting and making sure that 0 does not appear twice    
    if isempty(find(qc==0.5))==0    % then q=0.5 is a cut-off
        qc1 = [qc2(1:(end+1)/2);qc2((end+1)/2:end)];
    else
        qc1 = qc2;
    end
    q = qc1(end/2+1:end);
    
    k = length(q);
    xi = nines(T,k);    % variable used in constructing test statistic
    rhodiff = nines(k,1);
    for jj = 1:k;
        temp1 = find( (X<=quantile(X,1-qc(jj))).*(Y<=quantile(Y,1-qc(jj))) );
        temp2 = find( (X>=quantile(X,qc(jj))).*(Y>=quantile(Y,qc(jj))) );
        out1(jj,1) = corrcoef12(X(temp1),Y(temp1));
        out1(end+1-jj,1) = corrcoef12(X(temp2),Y(temp2));
        rhodiff(jj) = corrcoef12(X(temp2),Y(temp2)) - corrcoef12(X(temp1),Y(temp1));      % computing the test stat
        
        x1plus = (X-mean(X(temp2)))./std(X(temp2));
        x2plus = (Y-mean(Y(temp2)))./std(Y(temp2));
        x1minus = (X-mean(X(temp1)))./std(X(temp1));
        x2minus = (Y-mean(Y(temp1)))./std(Y(temp1));
        
        xi(:,jj) = T/size(temp2,1)*(x1plus.*x2plus - corrcoef12(X(temp2),Y(temp2))).*( (X>=quantile(X,qc(jj))).*(Y>=quantile(Y,qc(jj))) ) - ...
            T/size(temp1,1)*(x1minus.*x2minus - corrcoef12(X(temp1),Y(temp1))).*( (X<=quantile(X,1-qc(jj))).*(Y<=quantile(Y,1-qc(jj))) );
    end
    omegahat = newey_west(xi);      % computing the covariance matrix
    teststat = T*rhodiff'*inv(omegahat)*rhodiff;
    pval = 1-chi2cdf(teststat,k);
end
if indic==0
    X = (X-mean(X))./std(X);
    Y = (Y-mean(Y))./std(Y);
    
    % firstly: ripping through to make sure I have enough data in each tail
    qc2 = qc;
    for jj=1:length(qc);
        temp1 = find( (X<=-qc(jj)).*(Y<=-qc(jj)) );
        temp2 = find( (X>=qc(jj)).*(Y>=qc(jj)) );
        if ((length(temp1)<3)+(length(temp2)<3))>0 % then not enough data in at least one tail
            qc2 = setdiff(qc2,qc(jj));
        end
    end

    % now I am sure that each cut-off has enough data
    qc2 = sortrows(unique([qc2;-qc2]));   % sorting and making sure that 0 does not appear twice    
    if isempty(find(qc==0))==0    % then c=0 is a cut-off
        qc1 = [qc2(1:(end+1)/2);qc2((end+1)/2:end)];
    else
        qc1 = qc2;
    end
    c = qc1(end/2+1:end);
    
    k = length(c);
    xi = nines(T,k);    % variable used in constructing test statistic
    out1 = nines(2*k,1);		
    rhodiff = nines(k,1);
    for jj = 1:k;
        temp1 = find( (X<=-c(jj)).*(Y<=-c(jj)) );
        out1(jj,1) = corrcoef12(X(temp1),Y(temp1));
        temp2 = find( (X>=c(jj)).*(Y>=c(jj)) );      
        out1(end+1-jj,1) = corrcoef12(X(temp2),Y(temp2));
        rhodiff(jj) = corrcoef12(X(temp2),Y(temp2)) - corrcoef12(X(temp1),Y(temp1));      % computing the test stat
        
        x1plus = (X-mean(X(temp2)))./std(X(temp2));
        x2plus = (Y-mean(Y(temp2)))./std(Y(temp2));
        x1minus = (X-mean(X(temp1)))./std(X(temp1));
        x2minus = (Y-mean(Y(temp1)))./std(Y(temp1));
        
        xi(:,jj) = T/size(temp2,1)*(x1plus.*x2plus - corrcoef12(X(temp2),Y(temp2))).*((X>=c(jj)).*(Y>=c(jj))) - ...
            T/size(temp1,1)*(x1minus.*x2minus - corrcoef12(X(temp1),Y(temp1))).*((X<=-c(jj)).*(Y<=-c(jj)));
    end
    omegahat = newey_west(xi);      % computing the covariance matrix
    teststat = T*rhodiff'*inv(omegahat)*rhodiff;
    pval = 1-chi2cdf(teststat,k);
end
