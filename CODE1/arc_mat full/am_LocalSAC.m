classdef am_LocalSAC<am_LocalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display local estimates of SAC model
% ----------------------------------------------------------------------
% see also: LocalModel, GlobalSDM, GlobalSEM, GlobalSAR
% ----------------------------------------------------------------------
   properties
   end

   methods

%   Name: LocalSAC
%   Function: Construction function of class 'LocalSAC'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a LocalSAC object 

        function obj=am_LocalSAC(dataObj,y,x)
           obj = obj@am_LocalModel(dataObj);
           obj.ModelName = 'Local SAC Estimates';
           obj.HFig = 0;
           obj.HEnableSyn = true;
        end

%   Name: callModel
%   Function: Call underlying sac function to compute SAC estimates
%   Input: obj = current LocalSAC object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SAC estimates

       function sresults = callModel(obj,y,x,W)
           sresults = sac(y,x,W,W);
       end
       
%   Name: updatePlot
%   Function: Since local models are always computed after user selection of
%   map polygons, we do not present LocalSAC result with 'onPlot' function, and use 
%   'updatePlot' to present LocalSAC result when selections are made. 
%   Input: obj = current LocalSAC object
%   Output; Produce a tabular view of local SAC estimates

       function updatePlot(obj)
           updatePlot@am_LocalModel(obj);
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
