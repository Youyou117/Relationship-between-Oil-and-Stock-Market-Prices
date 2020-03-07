function [LLF, h, likelihoods] = fattailed_garchlikelihood2(parameters , data , p , q, errortype)

[r,c]=size(parameters);
if c>r
   parameters=parameters';
end


parameters(find(parameters <= 0)) = realmin;

constp=parameters(1);
archp=parameters(2:p+1);
garchp=parameters(p+2:p+q+1);
if errortype==2 | errortype==3
   nu = parameters(p+q+2);
   parameters = parameters(1:p+q+1);
elseif errortype==4
   nu = parameters(p+q+2:p+q+3);
   parameters = parameters(1:p+q+1);
end


if isempty(q)
   m=p;
else
   m  =  max(p,q);   
end

Tau=size(data,1);
stdEstimate =  std(data,1);                      
data        =  [stdEstimate(ones(m,1)) ; data];  
T           =  size(data,1);                    


h  =  data(1).^2;
h  =  [h(ones(m,1)) ; zeros(T-m,1)];   % Pre-allocate the h(t) vector.

h=garchcore(h,data,parameters,p,q,m,T);
%for t = (m + 1):T
%   h(t) = parameters' * [1 ; data(t-(1:p)).^2;  h(t-(1:q)) ];
%end



Tau = T-m;
LLF = 0;
t = (m + 1):T;
if errortype == 1
   LLF  =  sum(log(h(t))) + sum((data(t).^2)./h(t));
   LLF  =  0.5 * (LLF  +  (T - m)*log(2*pi));
elseif errortype == 2
   LLF = Tau*gammaln(0.5*(nu+1)) - Tau*gammaln(nu/2) - Tau/2*log(pi*(nu-2));
   LLF = LLF - 0.5*sum(log(h(t))) - ((nu+1)/2)*sum(log(1 + (data(t).^2)./(h(t)*(nu-2)) ));
   LLF = -LLF;
elseif errortype == 3
   Beta = (2^(-2/nu) * gamma(1/nu)/gamma(3/nu))^(0.5);
   LLF = (Tau * log(nu)) - (Tau*log(Beta)) - (Tau*gammaln(1/nu)) - Tau*(1+1/nu)*log(2);
   LLF = LLF - 0.5 * sum(log(h(t))) - 0.5 * sum((abs(data(t)./(sqrt(h(t))*Beta))).^nu);
   LLF = -LLF;
else 
    x = data(t)./sqrt(h(t));
    c = gamma((nu(1)+1)/2)./(sqrt(pi*(nu(1)-2)).*gamma(nu(1)/2));
    a = 4*nu(2).*c.*((nu(1)-2)./(nu(1)-1));
    b = sqrt(1 + 3*nu(2).^2 - a.^2);
    
    logc = gammaln((nu(1)+1)/2) - gammaln(nu(1)/2) - 0.5*log(pi*(nu(1)-2));
    logb = 0.5*log(1 + 3*nu(2).^2 - a.^2);
    
    find1 = (x<(-a./b));
    find2 = (x>=(-a./b));
    LL1   = logb + logc - (nu(1)+1)/2.*log(1+1./(nu(1)-2).*((b.*x+a)./(1-nu(2))).^2) - 0.5*log(h(t));
    LL2   = logb + logc - (nu(1)+1)/2.*log(1+1./(nu(1)-2).*((b.*x+a)./(1+nu(2))).^2) - 0.5*log(h(t));
    LLF   = sum(LL1(find1)) + sum(LL2(find2));
    LLF   = -LLF;    
end


if nargout > 2
   likelihoods=zeros(size(T));
   if errortype == 1
      likelihoods = 0.5 * ((log(h(t))) + ((data(t).^2)./h(t)) + log(2*pi));
      likelihoods = -likelihoods;
   elseif errortype == 2
      likelihoods = gammaln(0.5*(nu+1)) - gammaln(nu/2) - 1/2*log(pi*(nu-2))...
         - 0.5*(log(h(t))) - ((nu+1)/2)*(log(1 + (data(t).^2)./(h(t)*(nu-2)) ));
      likelihoods = -likelihoods;
  elseif errortype == 3
      Beta = (2^(-2/nu) * gamma(1/nu)/gamma(3/nu))^(0.5);
      likelihoods = (log(nu)/(Beta*(2^(1+1/nu))*gamma(1/nu))) - 0.5 * (log(h(t))) ...
         - 0.5 * ((abs(data(t)./(sqrt(h(t))*Beta))).^nu);
      likelihoods = -likelihoods;
  elseif errortype == 4 
    likelihoods(find1) = LL1(find1);
    likelihoods(find2) = LL2(find2);
    likelihoods = -likelihoods';
  end
end

h=h(t);
