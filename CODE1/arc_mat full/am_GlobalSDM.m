classdef am_GlobalSDM<am_GlobalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display global estimates of Spatial Durbin Model (SDM)
% ----------------------------------------------------------------------
% see also: GlobalModel, GlobalSEM, GlobalSAC, GlobalSAR
% ----------------------------------------------------------------------

   properties
   end

   methods
       
%   Name: GlobalSDM
%   Function: Construction function of class 'GlobalSDM'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a GlobalSDM object 

        function obj =  am_GlobalSDM(dataObj,y,x)
           obj = obj@am_GlobalModel(dataObj);
           obj.ModelName = 'Global SDM Estimates';
           obj.HFig = 0;
           obj.onPlot;           
        end

%   Name: callModel
%   Function: Call underlying Sdm function to compute
%   spatial regression model estimates
%   Input: obj = current GlobalSDM object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SDM model estimates
          
       function sresults = callModel(obj,y,x,W)
           sresults = sdm(y,x,W);
       end      
       
%   Name: onPlot
%   Function: Produce a tabular view of SDM model estimates
%   Input: obj = current GlobalSDM object
%   Output; 
       
       function onPlot(obj)
            onPlot@am_GlobalModel(obj);
       end
       function estimates = getEstimate(obj,y,x,W,sre)
            for q = 1:size(sre.beta,1)
                beta(q) = roundoff(sre.beta(q),4);
            end
              rho = roundoff(sre.rho,4);
            [n junk] = size(W);
            IN = eye(n); 
            xmat = [x W*x(:,2:end)];
            estimates = y-(IN - rho*W)\(xmat*beta');
       end

   end
end 
