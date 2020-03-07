classdef am_LocalSAR <am_LocalModel
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Display local estimates of SAR model
% ----------------------------------------------------------------------
% see also: LocalModel, GlobalSDM, GlobalSEM, GlobalSAC
% ----------------------------------------------------------------------

   properties
   end

   methods

%   Name: LocalSAR
%   Function: Construction function of class 'LocalSAR'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a LocalSAR object 

       function obj=am_LocalSAR(dataObj,y,x)
           obj = obj@am_LocalModel(dataObj);
           obj.ModelName = 'Local SAR Estimates';
           obj.HFig = 0;
           obj.HEnableSyn = true;
       end
       
       
%   Name: callModel
%   Function: Call underlying sar function to compute SAR estimates
%   Input: obj = current LocalSAR object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the SAR estimates

       function sresults = callModel(obj,y,x,W)
           sresults = sar(y,x,W);
       end
       
%   Name: updatePlot
%   Function: Since local models are always computed after user selection of
%   map polygons, we do not present LocalSAR result with 'onPlot' function, and use 
%   'updatePlot' to present LocalSAR result when selections are made. 
%   Input: obj = current LocalSAR object
%   Output; Produce a tabular view of local SAR estimates

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
           estimates = y - (IN-rho*W)\(x*beta');
       end

   end
  end

