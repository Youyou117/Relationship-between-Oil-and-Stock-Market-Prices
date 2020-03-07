classdef am_BaseMap <am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Produce a map or a choropleth map with user-specified variables          
% ----------------------------------------------------------------------
% see also: Representation
% ----------------------------------------------------------------------


   properties
       poly = [];
       obs_selected = [];  
       new_variable = []; 
       new_missing = []; 
       map_color = [];
   end
%   poly = a structure variable with handles to the map polygons, default = []
%               poly(i).handles(k) = a handle to each of (nobs=npoly) polygon regions and its k-parts
%                                   (these are patch objects)
%               poly(1).fig_handle = a handle to a figure containing the map
%                                   (use: set(poly(1).fig_handle,'Visible','on') to see the map)
%   obs_selected = a vector of 0's and 1's for selected and non-selected observations
%               0 = selected, 1 = non-selected, default = [] 
%   new_variable = a temporay vector of observation values, default = []
%   new_missing = a temporary vector of 0's and 1's for missing and non-missing observations
%               0 = missing, 1 = non-missing, default = []
%               (produces white map polygons for missing
%               values)

   methods
       
%   Name: BaseMap
%   Function: Construction function of class 'BaseMap'
%   Input: dataObj = a DataSource object to be mapped
%   Output: obj = a BaseMap object

       function obj = am_BaseMap(dataObj)

           obj = obj@am_Representation(dataObj); %call construction function of superior class 'Representation'
           obj.onPlot; % plot map
           
       end
       
%   Name: onPlot
%   Function: Produce a map 
%   Input: obj = current BaseMap object 
%   Output: poly is assigned with a structure variable returned by
%   read_shape() and contains data matrix to be mapped
      
       function onPlot(obj)
            
            results = obj.DataObj.results; %trade off between time and space; '.' operator is time consuming; therefore, we simply use local variables here
            poly  = make_map(results); %call underlying make_map function, returning a structure variable with handles to the map polygons
            
            hd = uicontextmenu; %create context menu
            obj.HEnableCm = uimenu(hd,...
            'Label','Synchronization',...
            'Checked','on',...
            'Callback', @(src,evnt)enableSyn(obj,src,evnt)); %create context menu item 'allow sychronization on this view'
            obj.HDisableCm = uimenu(hd,...
            'Label','Non-synchronization',...
            'Checked','off',...
            'Callback', @(src,evnt)disableSyn(obj,src,evnt)); %create context menu item 'do not allow sychronization on this view'
           set(obj.HFig,'uicontextmenu',hd); 
           
           vflag = results.vflag; 
           nvarsm = results.nvarsm;
           vnames = results.vnames;
           mname = ['Map of ' results.vnames(1,:)];
           if results.mapmenu == 0 %show menu toolbar
                set(poly(1).fig_handle, ...
                 'Position',[10 100 800 800], ... % [left bottom width height]
                 'MenuBar','none', ...
                 'NumberTitle','off', ...
                 'Name',mname,...
                 'Toolbar','figure');
            elseif results.mapmenu == 1
                set(poly(1).fig_handle, ...
                 'Position',[10 100 800 800], ...
                 'NumberTitle','off', ...
                 'Name',mname,...
                 'Toolbar','figure');
           end

            figure(poly(1).fig_handle);
            axis equal;
 
            if ~isempty(results.mapcolors) %if results.mapcolors is not empty, we are creating a colored map
                hold on; %hold on the figure, otherwise the whole map will be drawn everytime we add a polygon onto the screen 
                for i=1:results.npoly;
                    for k=1:results.nparts(i);
                        set(poly(i).handles(k),'FaceColor', results.mapcolors(i,:),'Visible','on'); %set polygon color
                    end;
                end;
                hold off;
            else
                hold on;
                thandles = zeros(results.npoly,1);
                for i=1:results.npoly;
                    for k=1:results.nparts(i);
                        set(poly(i).handles(k),'FaceColor',[1 1 1],'Visible','on'); %draw white polygons
                    end;
                end;
                hold off;
            end

        spop = uicontrol('Style', 'popup',...
            'String', 'Selections|select regions|add points to regions| subtract points from regions|de-select regions',...
            'Position', [10 5 200 20],...
             'Callback',@(src,evnt)mapselection(obj,src,evnt)); %map selection popup-menu

        qpop = uicontrol('Style', 'popup',...
            'String', 'Quit|Close/Exit',...
            'Position', [551 5 75 20],...
            'Callback',@(src,evnt)quitselection(obj,src,evnt)); %quit program popup-menu


        vlist = 'Variable|'; %create a list of variables
        if vflag == 1
            for i=1:nvarsm;
                vlist = [vlist vnames(i,:) '|'];
            end;
        elseif vflag == 0
            for i=1:nvarsm;
                vlist = [vlist vnames(i+1,:) '|'];
            end;
        end;

       vpop = uicontrol('Style', 'popup',...
            'String', vlist,...
            'Position', [331 5 150 20],...
            'Callback', @(src,evnt)variableselection(obj,src,evnt)); %variable selection popup-menu
        
       set(poly(1).fig_handle,'Visible','on'); %set all polygons to be visible
       
       obj.poly  = poly; %assign local variables' value to object properties
       obj.HFig = poly(1).fig_handle;
       obj.DataObj.results = results;
       
       end
       
%   Name: updatePlot
%   Function: Update current map according to user operations
%   Input: obj = current BaseMap object 
%   Output: 

       function updatePlot(obj)
           
           results = obj.DataObj.results; %trade off between time and space; '.' operator is time consuming; therefore, we simply use local variables here
           npoly = results.npoly;
           nparts = results.nparts;
           poly = obj.poly;
           
           figure(obj.HFig); %make the BaseMap object 'current figure'
           
           if ~isempty(results.mapcolors) %determine map polygons' color
                hold on;
                for i=1:npoly
                    for k=1:nparts(i)
                        set(poly(i).handles(k),'FaceColor',results.mapcolors(i,:),'Visible','on');
                    end;
                end;
           else
                hold on;
                thandles = zeros(npoly,1);
                obj.DataObj.results.mapcolors = ones(npoly,3);
                for i=1:npoly;
                    for k=1:nparts(i)
                        set(poly(i).handles(k),'FaceColor',[1 1 1],'Visible','on');
                    end;
                end;
           end
           
           obs_selected = results.obs_selected;
           obs_selected1 = results.obs_selected1; % previously selected map polygons
           obs_selected = results.obs_selected; %currently selected map polygons
           
           for i=1:length(obs_selected1) %redraw the boundary of previously selected map polygons
                p = obs_selected1(i);
                for k = 1:nparts(p)
                    set(poly(p).handles(k),'EdgeColor',[0 0 0],'LineWidth',1);
                end;
           end;
    
           for i=1:length(obs_selected) 
               p = obs_selected(i);
               for k = 1:nparts(p)
                    set(poly(p).handles(k),'EdgeColor',[0 0 0],'LineWidth',2); %draw the boundary of currently selected polygons with thick lines
               end;
           end;
          
           highlightpoly1 =obj.DataObj.results.highlightpoly1; %highlighted map polygons correspond to datapoints selected in the linked histogram
           if ~isempty(highlightpoly1)
                for i = 1:length(highlightpoly1)
                    p = highlightpoly1(i);
                    for k = 1:nparts(p)
                        set(poly(p).handles(k),'FaceColor',results.mapcolors(p,:));
                    end;
                end
           end
           if ~isempty(obj.DataObj.results.highlightpoly)
                highlightpoly = obj.DataObj.results.highlightpoly;
                for i = 1:length(highlightpoly)
                    p = highlightpoly(i);
                    for k = 1:nparts(p)
                        set(poly(p).handles(k),'FaceColor',[0.5 0.5 0.5]);
                    end;
                end
           end
           hold off;
           
       end
%   Name: mapselection
%   Function: Handle selection of map polygons
%   Input: obj = current BaseMap object
%           src = a source object 
%           evnt = a event object
%   Output: selection information in the DataSource object is updated     

       function mapselection(obj,src,evnt)
        
            temp=zoom;
            tst = get(temp,'Enable'); %set off 'zoom in' mode
            if strcmp(tst,'on')
                zoom;
            end

            temp = pan;
                tst = get(temp,'Enable'); %set off 'pan' mode
            if strcmp(tst,'on')
                pan;
            end;
            
            val = get(src,'Value'); %get current selection mode
            figure(obj.poly(1).fig_handle);
            variable = obj.DataObj.results.variable; %avoid using '.'operator in future processing
            vindex = obj.DataObj.results.vindex;
            missing = obj.DataObj.results.missing;
            obj.DataObj.results.highlightpoly1 = [];
            obj.DataObj.results.highlightpoly = [];
            obj.DataObj.results.highlightclass = 0;
    
            if val ==2 %select by rectangle
                if obj.DataObj.results.Hwarning ~= 0
                    set(obj.DataObj.results.Hwarning,'Visible','off');
                end;
                k = waitforbuttonpress;
                point1 = get(gca,'CurrentPoint');    % button down detected
                finalRect = rbbox;                   % return figure units
                point2 = get(gca,'CurrentPoint');    % button up detected
                point1 = point1(1,1:2);              % extract x and y
                point2 = point2(1,1:2);
                p1 = min(point1,point2);             % calculate locations
                offset = abs(point1-point2);         % and dimensions
                x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
                y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];

                xlim1 = min(x);
                xlim2 = max(x);
                ylim1 = min(y);
                ylim2 = max(y);

                maxx = xlim2;
                maxy = ylim2;
                minx = xlim1;
                miny = ylim1;

                figure(obj.poly(1).fig_handle);

                xmin = obj.DataObj.results.xmin;
                xmax = obj.DataObj.results.xmax;
                ymin = obj.DataObj.results.ymin;
                ymax = obj.DataObj.results.ymax;
                
                obs_selectedt = [];
                npoly = obj.DataObj.results.npoly;
                for i=1:npoly
                    ind = find(xmin(i) > xlim1 & xmax(i) < xlim2 & ymin(i) > ylim1 & ymax(i) < ylim2, 1);
                    if isempty(ind)
                        %       for k = 1:obj.DataObj.results.nparts(i);
                        %       set(obj.poly(i).handles(k),'FaceColor',[1 1 1]); % turns off previous selections
                        %       end;
                    else
                        obs_selectedt = [obs_selectedt
                        i];
                    end;
               end;
               obj.obs_selected = obs_selectedt;
               obj.new_variable =  variable(obj.obs_selected,:); 
               obj.new_missing = missing(obj.obs_selected,1);

               if isempty(obj.obs_selected); % the user botched a selection
                    obj.obs_selected = [];
                    obj.new_variable = variable;
                    obj.new_missing = missing;
               end;

               obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected;
               obj.DataObj.results.obs_selected = obj.obs_selected; 

               obj.DataObj.results.cvariable = obj.new_variable(:,vindex); % we send down the current variable to create the legend
               obj.DataObj.results.svariable = obj.new_variable; % contains all variables zoomed
               obj.DataObj.results.cmissing = obj.new_missing;

            elseif val ==3%select polygon by pointing

                if length(obj.DataObj.results.obs_selected) == obj.DataObj.results.npoly %all polygons have been selected
                    figure(obj.poly(1).fig_handle);
          
                    if obj.DataObj.results.Hwarning == 0
                        Hwarning=uicontrol('Style','text','Units','Normalized','Position',[0.1,0.1,1,0.035],'Backgroundcolor',[1 1 1], ...
                                'String',['You must first select a region'],'enable','inactive','FontSize',6,'HorizontalAlignment','left');
                        obj.DataObj.results.Hwarning = Hwarning;
                    else
                        set(obj.DataObj.results.Hwarning,'Visible','on');
                    end;
                    return;
                end;

            [xselect,yselect] = ginput;
            nselect = length(xselect);
            xc =obj.DataObj.results.xc;
            yc = obj.DataObj.results.yc;
            for i=1:nselect
                xi = xselect(i);
                yi = yselect(i);

                % find distance to latt-long centroids
                dist = (xi - xc).*(xi - xc) + (yi -yc).*(yi - yc);
                [mdist,ind] = min(dist);
                temp = find(obj.obs_selected == ind, 1);
                if isempty(temp)
                    obj.obs_selected = [obj.obs_selected
                                                ind];
                end
                obj.new_variable =  variable(obj.obs_selected,:); 
                obj.new_missing = missing(obj.obs_selected,1);
            end;
            obj.DataObj.results.cvariable = obj.new_variable(:,vindex); % we send down the current variable to create the legend
            obj.DataObj.results.svariable = obj.new_variable; % contains all variables zoomed
            obj.DataObj.results.cmissing = obj.new_missing;
            obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected;
            obj.DataObj.results.obs_selected = obj.obs_selected;   
    
            elseif val ==4 %delete single polygon from current selection 

            if length(obj.DataObj.results.obs_selected) == obj.DataObj.results.npoly
                figure(obj.poly(1).fig_handle);
                if obj.DataObj.results.Hwarning == 0
                    Hwarning=uicontrol('Style','text','Units','Normalized','Position',[0.1,0.1,1,0.035],'Backgroundcolor',[1 1 1], ...
                            'String',['You must first select a region'],'enable','inactive','FontSize',6,'HorizontalAlignment','left');
                    obj.DataObj.results.Hwarning = Hwarning;
                else
                    set(obj.DataObj.results.Hwarning,'Visible','on');
                end;
                return;
            end;

            [xselect,yselect] = ginput;
            nselect = length(xselect);

            pitch_obs = [];
            xc = obj.DataObj.results.xc;
            yc = obj.DataObj.results.yc;
            for i=1:nselect
                xi = xselect(i);
                yi = yselect(i);
                % find distance to latt-long centroids
                dist = (xi - xc).*(xi - xc) + (yi - yc).*(yi - yc);
                [mdist,ind] = min(dist);
                temp = find(obj.obs_selected == ind, 1);
                if ~isempty(temp)
                    pitch_obs = [pitch_obs
                                    ind];
                else
                    return;
                end    
            end;

            newobs = [];
            newvariable =[];
            newmissing = [];
            for i=1:length(obj.obs_selected)
                obsi = obj.obs_selected(i);
                ind = find(obsi == pitch_obs, 1);
                if isempty(ind)
                    newobs = [newobs
                                obsi];
                    newvariable = [newvariable
                                variable(obsi,:)];
                    newmissing = [newmissing
                                missing(obsi,1)];               
                else
                % do nothing
                end;
            end;
        obj.obs_selected = newobs;
        obj.new_variable = newvariable;
        obj.new_missing = newmissing;
        obj.DataObj.results.cvariable = obj.new_variable(:,vindex); % we send down the current variable to create the legend
        obj.DataObj.results.svariable = obj.new_variable; % contains all variables zoomed
        obj.DataObj.results.cmissing = obj.new_missing;
        obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected;
        obj.DataObj.results.obs_selected = obj.obs_selected; 
    elseif val==5%botch all selection
        if obj.DataObj.results.Hwarning ~= 0
            set(obj.DataObj.results.Hwarning,'Visible','off');
        end;
        obj.obs_selected = [];
        obj.new_variable = [];
        obj.new_missing = [];
        obj.DataObj.results.cvariable = variable(:,vindex); % we send down the current variable to create the legend
        obj.DataObj.results.svariable = variable; % contains all variables zoomed
        obj.DataObj.results.cmissing = missing;
        obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected;
        obj.DataObj.results.obs_selected = obj.obs_selected;  

        figure(obj.HFig);
    end   
    obj.update%call update function
    
       end   
       
%   Name: variableselection
%   Function: Handle selection of variable under invetstigation
%   Input: obj = current BaseMap object
%           src = a source object 
%           evnt = a event object
%   Output: selection information in the DataSource object is updated   

       function variableselection(obj,src,evnt)
          
           temp=zoom;
           tst = get(temp,'Enable'); %set off 'zoom in' mode
           if strcmp(tst,'on')
                zoom;
           end

            temp = pan;
            tst = get(temp,'Enable'); %set off 'pan' mode
            if strcmp(tst,'on')
                pan;
            end;
            
            obj.DataObj.results.highlightpoly1 = [];
            obj.DataObj.results.highlightpoly = [];
            obj.DataObj.results.highlightclass = 0;

            vnames = obj.DataObj.results.vnames;
            variable = obj.DataObj.results.svariable; % all variables possibly zoomed

            vindex = get(src,'Value');
            if vindex ~=1
                vindex = vindex - 1;
                vflag = obj.DataObj.results.vflag;
            end

            obj.DataObj.results.cvariable = variable(:,vindex);% current variable selection
            obj.DataObj.results.svariable = variable;
            obj.DataObj.results.vindex = vindex;
            % uses results.cvariable

            if vflag == 1
                set(obj.poly(1).fig_handle,'Name',['Map of ' obj.DataObj.results.vnames(vindex,:)]);
            elseif vflag == 0
                set(obj.poly(1).fig_handle,'Name',['Map of ' obj.DataObj.results.vnames(vindex+1,:)]);
            end;
            
            obj.update
            
       end
       
%   Name: quitselection
%   Function: Quit the program
%   Input: obj = current BaseMap object
%           src = a source object 
%           evnt = a event object
%   Output:    

       function quitselection(obj,src,evnt)
              close all;
       end
       
%   Name: update
%   Function: Send out 'map2graph' event
%   Input: obj = current BaseMap object
%   Output: a 'map2graph' event is sent to produce dynamically linking and brushing 

       function update(obj)
           notify(obj,'map2graph');%send out  'map2graph' event
       end
   end
end 
