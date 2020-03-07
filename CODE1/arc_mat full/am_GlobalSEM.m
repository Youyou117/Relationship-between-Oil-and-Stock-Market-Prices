classdef am_GlobalSEM<am_GlobalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display global estimates of Spatial Error Model (SEM)
% ----------------------------------------------------------------------
% see also: GlobalModel, GlobalSDM, GlobalSAC, GlobalSAR
% ----------------------------------------------------------------------

   properties
   end

   methods
       
%   Name: GlobalSEM
%   Function: Construction function of class 'GlobalSEM'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a GlobalSEM object 

       function obj =  am_GlobalSEM(dataObj,y,x)
           obj = obj@am_GlobalModel(dataObj);
           obj.ModelName = 'Global SEM Estimates';
           obj.onPlot;
           
       end
       
%   Name: callModel
%   Function: Call underlying sdm function to compute SDM model estimates
%   Input: obj = current GlobalSEM object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SEM model estimates

       function sresults = callModel(obj,y,x,W)
           sresults = sem(y,x,W);
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
           estimates = y - (x*beta');
       end
   end
end 