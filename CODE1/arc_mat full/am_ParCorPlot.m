classdef am_ParCorPlot<am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Produce a Parallel Coordinate Plot with user-specified variables 
% ----------------------------------------------------------------------
% see also: Representation
% ----------------------------------------------------------------------

   properties
        LineHd;
        Highlight;
   end
%   LineHd = a vector stores the handlers of polylines
%   A point in n-dimensional space is represented as a polyline with 
%	vertices on the parallel axes; the position of the vertex on the  
%   i-th axis corresponds to the i-th coordinate of the point.


   methods
       
%   Name: ParCorPlot
%   Function: Construction function of class 'ParCorPlot'
%   Input: dataObj = a DataSource object to be mapped
%   Output: obj = a ParCorPlot object 

       function obj = am_ParCorPlot(dataObj)
           obj = obj@am_Representation(dataObj);
           obj.LineHd = [];
           obj.onPlot
           hd = uicontextmenu;
            obj.HEnableCm = uimenu(hd,...
            'Label','Synchronization',...
            'Checked','on',...
            'Callback',@(src,evnt)enableSyn(obj,src,evnt));
            obj.HDisableCm = uimenu(hd,...
            'Label','Non-synchronization',...
            'Checked','off',...
            'Callback',@(src,evnt)disableSyn(obj,src,evnt));
           set(obj.HFig,'uicontextmenu',hd);
       end
       
%   Name: onPlot
%   Function: Produce a parallel coordinate plot 
%   Input: obj = current ParCorPlot object 
%   Output: Produce a moran scatterplot       

       function onPlot(obj)
           cmap = obj.DataObj.results.cmap;
            nbc = obj.DataObj.results.nbc;

            variable = obj.DataObj.results.variable; % we access this in case the user has 
                              % changed to a subset of the sample data 
            missing = obj.DataObj.results.cmissing;   % missing observations in the face of a zoom

            vindex = obj.DataObj.results.vindex;
            vflag = obj.DataObj.results.vflag;
            vnames = obj.DataObj.results.vnames;
            if vflag == 1
                vname = vnames(vindex,:);
            elseif vflag == 0
                vname = vnames(vindex+1,:);
            end;

            svec = get(0,'ScreenSize');
            if svec(3) > 1300
                width = 400; height = 600;
            elseif svec(3) > 1000
                width = 450; height = 650;
            elseif svec(3) == 800
                error(' you need a higher screen resolution than 800x600 to use arc_map');
            end;
            if obj.HFig == 0
                if obj.DataObj.results.legendmenu == 0
                    obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width height], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Parallel Coordinate Plot', ...
                      'MenuBar','none', ...
                      'Toolbar','figure');
                elseif obj.DataObj.results.legendmenu == 1
                    obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width height], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Parallel Coordinate Plot',...
                      'Toolbar','figure');
                end;
            else
                figure(obj.HFig);
            end;
            hc = colormap(cmap);

            spop1 = uicontrol('Style', 'popup',...
                    'String', 'select|clear selection',...
                    'Position', [2 5 150 20],...
                    'Callback',@(src,evnt)pcpselection(obj,src,evnt));
                    
           max1 = [];
           min1 = [];
           scale = [];
           for i = 1:size(variable,2)
               max1(i) = max(variable(:,i));
               min1(i) = min(variable(:,i));
               scale(i) = max1(i)-min1(i);
               pos = [i 100];
               text('Position',pos,'String',obj.DataObj.results.vnames(i,:),'Rotation',30);
               k = 1;
               tempX = [];
               tempY = [];
               tempX(k) = i;
               tempY(k) = 0;
               k = k+1;
               tempX(k) = i;
               tempY(k) = 100;
               line('XData',tempX,'YData',tempY,'Color','b','LineWidth',3);
           end
           for i=1:obj.DataObj.results.npoly
               tempX = [];
               tempY = [];
               tempX = 1:size(variable,2);
               for j = 1:length(tempX)
                    temp = (variable(i,j)-min1(j))/scale(j)*100;
                    tempY = [tempY temp]; 
               end
               obj.LineHd(i) = line('XData',tempX,'YData',tempY,'Color','k');
           end
           xlabel('Dimensions');
           ylabel('Relative Position');
       end
       
%   Name: updatePlot
%   Function: Update current parallel coodrinate plot
%   Input: obj = current ParCorPlot object
%   Output: ParCorPlot object is updated to highlight the selected data
%   polyline

       function updatePlot(obj)
           figure(obj.HFig);
           hold on;
           obs_selected = obj.DataObj.results.obs_selected;
           obs_selected1 = obj.DataObj.results.obs_selected1;
           lines = obj.LineHd;
           for i=1:length(obs_selected1)
                p = obs_selected1(i);  
                set(lines(p),'Color',[0 0 0],'LineWidth',1);
           end
           for i=1:length(obs_selected)
                p = obs_selected(i);
                set(lines(p),'Color',[0 1 1],'LineWidth',3);
           end
            highlightpoly1 =obj.DataObj.results.highlightpoly1;
            highlightpoly = obj.DataObj.results.highlightpoly;
            if ~isempty(highlightpoly1)
                for i = 1:length(highlightpoly1)
                    p = highlightpoly1(i);
                    set(lines(p),'Color',[0 0 0],'LineWidth',1);
                end
            end
            if ~isempty(highlightpoly)
                for i = 1:length(highlightpoly)
                    p = highlightpoly(i);
                    set(lines(p),'Color',[1 1 0],'LineWidth',3);
                end
            end
            hold off;
       end
       
%   Name: pcpselection
%   Function: Handle rectangle selection on the parallel coordinate plot
%   Input: obj = current ParCorPlot object
%           src = a source object 
%           evnt = a event object
%   Output: The information about selected data points (obs_selected) in
%   the DataSource object is updated

       function pcpselection(obj,src,evnt)
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
        val = get(src,'Value');
        figure(obj.HFig);
        variable = obj.DataObj.results.variable;
        vindex = obj.DataObj.results.vindex;
        missing = obj.DataObj.results.missing;

        obj.DataObj.results.highlightpoly1 = [];
        obj.DataObj.results.highlightpoly = [];
        obj.DataObj.results.highlightclass = 0;
        if val ==1
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

            obs_selected = [];
            lines = obj.LineHd;
            for i = 1:length(lines)
                tempx = get(lines(i),'XData');
                tempy = get(lines(i),'YData');
                [xi yi] = curveintersect(x,y,tempx,tempy);
                if ~isempty(xi)
                    obs_selected = [obs_selected
                                          i];                    
                end
            end 
  
            new_variable = variable(obs_selected,:);
            new_missing = missing(obs_selected,1); 
    
            obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected; 
            obj.DataObj.results.obs_selected = obs_selected;
    elseif val ==2
            obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected; 
            obj.DataObj.results.obs_selected = [];
            new_variable = variable;
            new_missing = missing;
    end
    
    obj.DataObj.results.cvariable = new_variable(:,vindex); % we send down the current variable to create the legend
    obj.DataObj.results.svariable = new_variable; % contains all variables zoomed
    obj.DataObj.results.cmissing = new_missing;
    notify(obj,'graph2map');
       end
       
   end
end





      
