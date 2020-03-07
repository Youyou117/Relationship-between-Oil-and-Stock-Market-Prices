classdef am_Linkage<handle
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Linkage object handles the dynamically linking and brushing of
%       representation objects. Since Representation objects register
%       themselves in the Linkage object using a handler for every
%       Representation object, developers can access Representation objects
%       and perform operations on those objects. For example, developers
%       can call class functions and modify class properties using the 
%       Linkage class. In addition, the Linkage object detects
%       Representation objects' events and updates them in a pre-defined
%       sequence.	
% ----------------------------------------------------------------------
% see also: Representation 
% ----------------------------------------------------------------------

   properties
       mapObj;
       hLMap;
       mapNum;
       graphObj;
       hLGraph;
       graphNum;
       graph2Obj;
       hLGraph2;
       graphNum2;
   end
%   mapObj = a vector stores the handler of BaseMap object
%   hLMap = a vectore stores the handler of eventListener to operations on
%           BaseMap object
%   mapNum = a variable stores the number of maps in the application,
%       currently we only allow one map in each application. This variable is
%       reserved for future multi-map applications.
%   graphObj = a vector stores the handlers of Representation objects that
%           may determine the color theme in BaseMap object
%   hLGraph = a vectore stores the handler of eventListener to operations on
%           Representation objects that may determine the color theme of
%           BaseMap object
%   graphNum = a vectore stores the number of Representation objects that
%           may determine the color theme of BaseMap object
%   graph2Obj = a vector stores the handlers of Representation objects that
%           do not determine the color theme in BaseMap object
%   hLGraph2 = a vectore stores the handler of eventListener to operations on
%           Representation objects that do not determine the color theme of
%           BaseMap object
%   graphNum2 = a vectore stores the number of Representation objects that
%           do not determine the color theme of BaseMap object
   methods

%   Name: Linkage
%   Function: Construction function of 'Linkage' class
%   Input: Null
%   Output: obj = a Linkage object  
       function obj = am_Linkage
           obj.initialize;        
       end

%   Name: initialize
%   Function: Initialize local variables 
%   Input: obj = current Linkage object
%   Output: 

       function initialize(obj) %set default values
           obj.mapObj = []; 
           obj.hLMap = 0;
           obj.mapNum = 0;
           obj.graphObj = {};
           obj.graph2Obj = {};
           obj.hLGraph = {};
           obj.hLGraph2 = {};
           obj.graphNum = 0;
           obj.graphNum2 = 0;
       end
       
%   Name: addMap
%   Function: Register BaseMap object in the Linkage object
%   Input: obj = current Linkage object
%           mapObj = BaseMap object to be registered
%   Output: the Linkage propertity 'mapObj' and 'hLMap' is updated

       function addMap(obj,mapObj)
           obj.mapObj = mapObj;
           obj.mapNum =1;
           obj.hLMap = addlistener(obj.mapObj,'map2graph',... %add listner to 'map2graph' event
            @(src,evnt)listenMap2Graph(obj,src,evnt));
       end

%   Name: addGraph
%   Function: Register Representation object that may determine BaseMap
%           objects' color theme in the Linkage object
%   Input: obj = current Linkage object
%           graph = Representation object to be registered
%   Output: the Linkage propertities: 'graphNum','graphObj', and 'hLGraph'
%   are updated

       function addGraph(obj,graph)
           obj.graphNum = obj.graphNum+1;
%            if obj.graphNum ~= 1 
%                 graph.HEnableSyn = false;
%                 set(graph.HEnableCm,'Checked','off') 
%                 set(graph.HDisableCm,'Checked','on')
%            end  
       
           %add listener to 'graph2map' event 
          obj.graphObj{obj.graphNum} =graph;
          obj.hLGraph{obj.graphNum} = {addlistener(obj.graphObj{obj.graphNum},'graph2map',... 
                                            @(src,evnt)listenGraph2Map(obj,src,evnt))};
          addlistener(obj.graphObj{obj.graphNum},'mapping',@(src,evnt)handleMapped(obj,src,evnt));
       end
      


%   Name: addGraph2
%   Function: Register Representation object that does not determine BaseMap
%           objects' color theme in the Linkage object
%   Input: obj = current Linkage object
%           graph = Representation object to be registered
%   Output: the Linkage propertities: 'graphNum2', 'graph2Obj', and
%   'hLGraph2' are updated

       function addGraph2(obj,graph)
          obj.graphNum2 = 1+ obj.graphNum2;           
          obj.graph2Obj{obj.graphNum2} =graph;
          %add listener to 'graph2map' event
          obj.hLGraph2{obj.graphNum2} = {addlistener(obj.graph2Obj{obj.graphNum2},'graph2map',... 
                                            @(src,evnt)listenGraph2Map(obj,src,evnt))};

       end
       
%   Name: listenMap2Graph
%   Function: Handle event 'map2graph' and brush all views
%   Input: obj = current Linkage object
%           src = a source object 
%           evnt = a event object
%   Output: All exisiting data views are refreshed

       function listenMap2Graph(obj,src,evnt)
           for i = 1:obj.graphNum
               if obj.graphObj{i}.HEnableSyn %update graphes that determine color theme of basemap
                    obj.graphObj{i}.updatePlot;
               end
           end
            if(obj.mapNum~=0)
                if (obj.mapObj.HEnableSyn) 
                    obj.mapObj.updatePlot;
                end
           end
           for i = 1:obj.graphNum2
               if obj.graph2Obj{i}.HEnableSyn
                    obj.graph2Obj{i}.updatePlot; %update other plots
               end
           end
       end

%   Name: listenGraph2Map
%   Function: Handle event 'graph2map' and brush all views
%   Input: obj = current Linkage object
%           src = a source object 
%           evnt = a event object
%   Output: All exisiting data views are refreshed 

       function listenGraph2Map(obj,src,evnt)
           for i = 1:obj.graphNum
               if obj.graphObj{i}.HEnableSyn
                    obj.graphObj{i}.updatePlot;
               end
           end
           if(obj.mapNum~=0)
                if (obj.mapObj.HEnableSyn) 
                    obj.mapObj.updatePlot;
                end
           end
           for i = 1:obj.graphNum2
               if obj.graph2Obj{i}.HEnableSyn
                    obj.graph2Obj{i}.updatePlot;
               end
           end
       end
       
 
   
        function handleMapped(obj,src,evnt)
            if ~isempty(obj.mapObj)
                if src.mapped ==1
                for i = 1:obj.graphNum
                    if (obj.graphObj{i} ~=src) && (obj.graphObj{i}.mapped == 1)
                       minVal = get(obj.graphObj{i}.chx,'Min'); 
                       set(obj.graphObj{i}.chx,'Value',minVal);
                   
                        obj.graphObj{i}.mapped=0;
                    end
  
                end
                src.DataObj.results.mapcolors = src.graph_color;
                         
                 if(obj.mapNum~=0)
                    if (obj.mapObj.HEnableSyn) 
                            obj.mapObj.updatePlot;
                    end
                 end
                else
                   src.DataObj.results.mapcolors =[];
                    if(obj.mapNum~=0)
                    if (obj.mapObj.HEnableSyn) 
                            obj.mapObj.updatePlot;
                    end
                 end
                end
            else
            end
        end
        end
   end
