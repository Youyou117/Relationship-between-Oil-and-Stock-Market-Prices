classdef am_LocalSDM<am_LocalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display local estimates of SDM model
% ----------------------------------------------------------------------
% see also: LocalModel, GlobalSAR, GlobalSEM, GlobalSAC
% ----------------------------------------------------------------------

   properties
   end

   methods
       
       
%   Name: LocalSDM
%   Function: Construction function of class 'LocalSDM'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a LocalSDM object 

         function obj=am_LocalSDM(dataObj,y,x)
           obj = obj@am_LocalModel(dataObj);
           obj.ModelName = 'Local SDM Estimates';
           obj.HFig = 0;
           obj.HEnableSyn = true;
         end
       
%   Name: callModel
%   Function: Call underlying sdm function to compute SDM estimates
%   Input: obj = current LocalSDM object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SDM estimates

       function sresults = callModel(obj,y,x,W)
           sresults = sdm(y,x,W);
       end
       
%   Name: updatePlot
%   Function: Since local models are always computed after user selection of
%   map polygons, we do not present LocalSDM result with 'onPlot' function, and use 
%   'updatePlot' to present LocalSDM result when selections are made. 
%   Input: obj = current LocalSDM object
%   Output; Produce a tabular view of local SDM estimates

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
            xmat = [x W*x(:,2:end)];
            estimates = y-(IN - rho*W)\(xmat*beta');
       end
   end
end 
