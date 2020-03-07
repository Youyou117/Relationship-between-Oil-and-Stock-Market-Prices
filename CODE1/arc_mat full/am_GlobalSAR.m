classdef am_GlobalSAR<am_GlobalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display global estimates of Spatial Autoregressive (SAR) model
% ----------------------------------------------------------------------
% see also: GlobalModel, GlobalSEM, GlobalSAC, GlobalSDM
% ----------------------------------------------------------------------

   properties
   end

   methods

%   Name: GlobalSAR
%   Function: Construction function of class 'GlobalSAR'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a GlobalSAR object 
       
       function obj =  am_GlobalSAR(dataObj)
           obj = obj@am_GlobalModel(dataObj);
           obj.ModelName = 'Global SAR Estimates';
          
           obj.HFig = 0;
           obj.onPlot;
           
       end

%   Name: callModel
%   Function: Call underlying SAR function to compute
%   spatial regression model estimates
%   Input: obj = current GlobalSAR object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SAR model estimates

       function sresults = callModel(obj,y,x,W)
           sresults = sar(y,x,W);
       end    
       
       function estimates = getEstimate(obj,y,x,W,sre)
            for q = 1:size(sre.beta,1)
                beta(q) = roundoff(sre.beta(q),4);
            end
            rho = roundoff(sre.rho,4);
            [n junk] = size(W);
            IN = eye(n); 
           estimates = y - (IN-rho*W)\(x*beta');
       end
       
%   Name: onPlot
%   Function: Produce a tabular view of SAR model estimates
%   Input: obj = current GlobalSAR object
%   Output; 

       function onPlot(obj)
            onPlot@am_GlobalModel(obj);
       end

   end
end 
