classdef am_GlobalSAC<am_GlobalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display global estimates of SAC model
% ----------------------------------------------------------------------
% see also: GlobalModel, GlobalSEM, GlobalSAR, GlobalSDM
% ----------------------------------------------------------------------


   properties
   end

   methods

%   Name: GlobalSAC
%   Function: Construction function of class 'GlobalSAC'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a GlobalSAC object 

        function obj =  am_GlobalSAC(dataObj)
           obj = obj@am_GlobalModel(dataObj);
           obj.ModelName = 'Global SAC Estimates';
           obj.HFig = 0;
           obj.onPlot;
           
        end
       
%   Name: callModel
%   Function: Call underlying SAC function to compute
%   spatial regression model estimates
%   Input: obj = current GlobalSAC object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SAC model estimates

       function sresults = callModel(obj,y,x,W)
           sresults = sac(y,x,W,W);
       end   
%   Name: onPlot
%   Function: Produce a tabular view of SAC model estimates
%   Input: obj = current GlobalSAC object
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
           estimates = y - (lN-W)\(x*beta');
       end
   end
end 
