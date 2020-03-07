classdef am_LocalSEM<am_LocalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display local estimates of SEM model
% ----------------------------------------------------------------------
% see also: LocalModel, GlobalSAR, GlobalSDM, GlobalSAC
% ----------------------------------------------------------------------

   properties
   end

   methods
       
%   Name: LocalSEM
%   Function: Construction function of class 'LocalSEM'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a LocalSEM object 

       function obj=am_LocalSEM(dataObj,y,x)
           obj = obj@am_LocalModel(dataObj);
           obj.ModelName = 'Local SEM Estimates';
           obj.HFig = 0;
           obj.HEnableSyn = true;
       end

%   Name: callModel
%   Function: Call underlying sem function to compute SEM estimates
%   Input: obj = current LocalSEM object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SEM estimates     
       
       
       function sresults = callModel(obj,y,x,W)
           sresults = sem(y,x,W);
       end

%   Name: updatePlot
%   Function: Since local models are always computed after user selection of
%   map polygons, we do not present LocalSEM result with 'onPlot' function, and use 
%   'updatePlot' to present LocalSEM result when selections are made. 
%   Input: obj = current LocalSEM object
%   Output; Produce a tabular view of local SEM estimates

       function updatePlot(obj)
           updatePlot@am_LocalModel(obj);
       end
       function estimates = getEstimate(obj,y,x,W,sre)
            for q = 1:size(sre.beta,1)
                beta(q) = roundoff(sre.beta(q),4);
            end
           estimates = y - (x*beta');
       end
   end
end 
