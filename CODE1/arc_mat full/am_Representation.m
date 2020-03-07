classdef am_Representation<handle
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
% The Representation classes are a set of classes that present 
% DataSource information in various ways.The set of Representation classes
% is organized in a hierarchical manner.All classes are inherited from the
% generic Repersentation class. The instantiation of Representation classes
% may be a choropleth map, statistical plot, or table showing modelling
% estimates.Different Representation objects share the same DataSource
% object to help researchers view data from different perspectives.The
% functions within Representation classes serve two general purposes:
% rendering the Graphical User Interface of the specific viewport, and
% handling users' operations on elements of the GUI.        
% ----------------------------------------------------------------------
% see also: BaseMap, Moranplot, PltDens, Histogram
% ----------------------------------------------------------------------
   properties
       DataObj
       HFig
       HEnableSyn
       HEnableCm
       HDisableCm
   end

%   DataObj = a vector stores the handler of DataSource object
%   HFig = a vector stores the handler of matlab Figure
%   HEnableSyn = a vector indicates whether the view synchronize with other
%                viewports, true = 'synchronization', false =
%                'non-synchronization'
%   HEnableCm = a vector indicates whether the context menu 'Enable Synchronization'
%                is checked or not, on = 'checked', off = 'unchecked'
%   HDisableCm = a vector indicates whether the context menu 'Disable Synchronization'
%                is checked or not, on = 'checked', off = 'unchecked'

   methods
       
%   Name: Representation
%   Function: Construction function of class 'Representation'
%   Input: dataObj = a DataSource object to be mapped
%   Output: obj = a Representation object 

       function obj = am_Representation(dataObj)
            obj.DataObj = dataObj;%store the handler of datasource object
            obj.HEnableSyn = true;
            obj.HFig = 0;
       end 
       
%   Name: onPlot
%   Function: Produce a view of the DataSource object 
%   Input: obj = current Representation object
%   Output: 

       function onPlot(obj)
       end

%   Name: updatePlot
%   Function: Update current view of the DataSource object
%   Input: obj = current Representation object
%   Output: 

       function updatePlot(obj)
       end

%   Name: enableSyn
%   Function: Enable synchronization among different views 
%   Input: obj = current Representation object
%           src = a source object 
%           evnt = a event object

       function enableSyn(obj,src,evnt)
           obj.HEnableSyn = true;
           set(obj.HEnableCm,'Checked','on')
           set(obj.HDisableCm,'Checked','off')
       end
       
%   Name: disableSyn
%   Function: Disable synchronization among different views 
%   Input: obj = current Representation object
%           src = a source object 
%           evnt = a event object

       function disableSyn(obj,src,evnt)
           obj.HEnableSyn = false;
           set(obj.HEnableCm,'Checked','off')
           set(obj.HDisableCm,'Checked','on')
       end
       
%   Name: initialize
%   Function: virtual function to be implemented by derived classes
       function initialize(obj)
       end
   end
   events
       map2graph
       graph2map
       mapping
% map2graph = a event representing operations on BaseMap and
%                synchronization among other views
% graph2map = a event representing operations on other views and
%               synchronization on the BaseMap
   end
end 
