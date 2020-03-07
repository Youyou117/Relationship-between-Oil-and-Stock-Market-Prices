classdef am_PltDens <am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Produce a density distribution curve for user-specified variable and
%       a density distribution curve for the selected subset. 
% ----------------------------------------------------------------------
% see also: Representation
% ----------------------------------------------------------------------
   properties
   end

   methods
       
%   Name: PltDens
%   Function: Construction function of class 'PltDens'
%   Input: dataObj = a DataSource object to be mapped
%   Output: obj = a PltDens object 

       function obj = am_PltDens(dataObj)
           obj = obj@am_Representation(dataObj);
           obj.HFig = 0;
           obj.onPlot;
       end
       
%   Name: onPlot
%   Function: Produce a density distribution plot 
%   Input: obj = current PltDens object
%   Output: Produce a general density distribution curve 
       function onPlot(obj)
           missing = obj.DataObj.results.cmissing;   % missing observations in the face of a zoom
           vindex = obj.DataObj.results.vindex;
            if isempty(obj.DataObj.results.obs_selected)
                variable = obj.DataObj.results.variable(:,vindex);
            elseif length(obj.DataObj.results.obs_selected)==1 %we may not create distribution plot with single data point
                variable = obj.DataObj.results.variable(:,vindex);
            else
            variable = obj.DataObj.results.variable(obj.DataObj.results.obs_selected,vindex); % we access this in case the user has 
                              % changed to a subset of the sample data 
            end
            vflag = obj.DataObj.results.vflag;
            vnames = obj.DataObj.results.vnames;
            if vflag == 1
                vname = vnames(vindex,:);
            elseif vflag == 0
                vname = vnames(vindex+1,:);
            end;


            svec = get(0,'ScreenSize');
            if svec(3) > 1300
                width = 400; height = 400;
            elseif svec(3) > 1000
                width = 450; height = 450;
            elseif svec(3) == 800
                error(' you need a higher screen resolution than 800x600 to use arc_map');
            end;
            
            if obj.HFig == 0
                if obj.DataObj.results.legendmenu == 0
                        obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width height], ... % [left bottom width height]
                        'NumberTitle','off', ...
                        'Name','Density Map', ...
                        'MenuBar','none', ...
                        'Toolbar','figure');
                elseif obj.DataObj.results.legendmenu == 1
                        obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width height], ... % [left bottom width height]
                        'NumberTitle','off', ...
                        'Name','Density Map',...
                        'Toolbar','figure');
                end;
                hd = uicontextmenu;
                obj.HEnableCm = uimenu(hd,...
                    'Label','Synchronization',...
                    'Checked','on',...
                    'Callback', @(src,evnt)enableSyn(obj,src,evnt));
                obj.HDisableCm = uimenu(hd,...
                    'Label','Non-synchronization',...
                    'Checked','off',...
                    'Callback', @(src,evnt)disableSyn(obj,src,evnt));

                set(obj.HFig,'uicontextmenu',hd);            
            else
                figure(obj.HFig);
            end;


            spop = uicontrol('Style', 'popup',...
                'String', 'Selections|select original density|deselect',...
                'Position', [0 5 200 20],...
                'Callback',@(src,evnt)densityselection(obj,src,evnt));

            % adjust the legend for missing values
            mindex = find(missing == 1);
            if isempty(mindex)
                wvariable = variable(mindex,1);
            else
                wvariable = variable;
            end;

            obj.DataObj.results.wvariable = wvariable; %local data of selected polygons
            mindex =find(obj.DataObj.results.globalcmissing ==1); %global missing represents the region with no data
            if isempty(mindex)
                variable1 = obj.DataObj.results.variable(mindex,vindex);
            else
                variable1 = obj.DataObj.results.variable(:,vindex);
            end
            obj.DataObj.results.variable1 = variable1; %gloabl data of all polygons
            

            [h1 f1 x1] = pltdens(wvariable);
            [h2 f2 x2] = pltdens(variable1);

            % force smoothness by doubling average of two
            % default bandwidths
            h = (h1+h2);

            [h1 f1 x1] = pltdens(wvariable,h);
            [h2 f2 x2] = pltdens(variable1,h);

            plot(x1,f1,x2,f2,'--','LineWidth',2);

            xlabel(vname);
       end

%   Name: updatePlot
%   Function: Update current moranplot
%   Input: obj = current PltDens object
%   Output: Produce one curve for the overall distribution of variable and
%   one curve for the distribution of selected subset

       function updatePlot(obj)
           obj.onPlot
       end

%   Name: densityselection
%   Function: Handle rectangle selection on the density plot
%   Input: obj = current PltDens object
%           src = a source object 
%           evnt = a event object
%   Output: The information about selected data points (obs_selected) in
%   the DataSource object is updated
       function densityselection(obj,src,evnt)
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
            t=cputime;
            val = get(src,'Value');
            figure(obj.HFig);
            vindex = obj.DataObj.results.vindex;
            variable = obj.DataObj.results.variable;            
            missing = obj.DataObj.results.missing;
            obj.DataObj.results.highlightpoly1 = [];
            obj.DataObj.results.highlightpoly = [];
            obj.DataObj.results.highlightclass = 0;
            obs_selected = [];

            if val ==2
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
                
                for i = 1:obj.DataObj.results.nobs
                    if  variable(i,vindex)<maxx && variable(i,vindex)>minx
                        obs_selected = [obs_selected
                                          i];                    
                    end
                end 
  
                new_variable = variable(obs_selected,:);
                new_missing = missing(obs_selected,1); 
    
                obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected; 
                obj.DataObj.results.obs_selected = obs_selected;
            else
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
