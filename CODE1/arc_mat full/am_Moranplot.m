classdef am_Moranplot<am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Produce a moranplot with user-specified variables and a color theme
%       for the coupled BaseMap object
% ----------------------------------------------------------------------
% see also: Representation
% ----------------------------------------------------------------------

   properties
       cvariable1; 
       WX;
       Patch;
       graph_color;
       chx;
   end
   properties (SetObservable)
       mapped;
   end
%   cvariable = a vector stores the local copy of selected variable values
%   WX = a vector stores the local copy of spatial contiguity matrix
%   Patch = a vector stores the handlers of scatter points

   methods
%   Name: Moranplot
%   Function: Construction function of class 'Moranplot'
%   Input: dataObj = a DataSource object to be mapped
%   Output: obj = a Moranplot object 
       function obj = am_Moranplot(dataObj)
           obj = obj@am_Representation(dataObj);
           latt = obj.DataObj.results.xc;
           long = obj.DataObj.results.yc;
           obj.DataObj.results.labels=0;     
           [j,W,j] = xy2cont(long,latt);%compute spatial weight matrix from lattitude and longtitude, users may define their own W matrix to relace this line
           %for example, users may use 'make_neighborsw(latt,long,n)'
           %function to create n-neighbors W matrix; 
           obj.DataObj.results.W =W;
           obj.initialize
           obj.HFig = 0;
           obj.graph_color = [];
           obj.mapped =0;
           obj.onPlot;
           obj.HEnableSyn=true;
       end
       
%   Name: initialize
%   Function: Initialize local variables and determine data points in four
%   quadrants 
%   Input: obj = current Moranplot object
%   Output: Quadrant information in the DataSource object is updated

       function initialize(obj)

           good = find(obj.DataObj.results.cmissing == 1);
           bad = find(obj.DataObj.results.cmissing == 0);
            [nobsm,nvarsm] = size(obj.DataObj.results.variable);
            tmp = mean(obj.DataObj.results.variable(good,:)); % adjust mean for missing values
            % the mean computed here is global mean
            svariable = matsub(obj.DataObj.variable,tmp); % center variables (transform)
            svariable(bad,:) = zeros(length(bad),nvarsm);
            obj.DataObj.results.svariable = svariable; % hold all variables (transformed)
            obj.DataObj.results.cvariable1 = svariable(:,obj.DataObj.results.vindex); % holds current variable selection (transformed)
            cvariable1 = obj.DataObj.results.cvariable1;
            obj.cvariable1 = svariable;
            WX = zeros(nobsm,nvarsm);
            WX =obj.DataObj.results.W(:,good)*svariable(good,:);

            % define the quadrants according to global means
            Q0=find(cvariable1 == 0 & WX(:,obj.DataObj.results.vindex)== 0);
            Q1=find(cvariable1>0 & WX(:,obj.DataObj.results.vindex)>0);
            Q2=find(WX(:,obj.DataObj.results.vindex)>0 & cvariable1<= 0);
            Q3=find(WX(:,obj.DataObj.results.vindex)<= 0 & cvariable1<= 0);
            Q4=find(cvariable1>0 & WX(:,obj.DataObj.results.vindex)<= 0);
            
            obj.DataObj.results.Q0 = Q0;
            obj.DataObj.results.Q1 = Q1;
            obj.DataObj.results.Q2 = Q2;
            obj.DataObj.results.Q3 = Q3;
            obj.DataObj.results.Q4 = Q4;
            obj.DataObj.results.WX = WX;
       end
       
%   Name: onPlot
%   Function: Produce a moranplot 
%   Input: obj = current Moranplot object 
%   Output: Produce a moran scatterplot

       function onPlot(obj)
           cmap = obj.DataObj.results.cmap;
           nbc = 4; %in moranplot, we have four 'categories', i.e., four quadrants
            vindex = obj.DataObj.results.vindex; %current selected variable index
            vflag = obj.DataObj.results.vflag;
            cvariable1 = obj.DataObj.results.cvariable1; 
            WX = obj.DataObj.results.WX;

            obs_selected = obj.DataObj.results.obs_selected; %current selected map polygons
            nobs = length(obs_selected);
            if nobs == 0
                nobs = obj.DataObj.results.npoly; 
                temp = 1:obj.DataObj.results.npoly;
                obs_selected = temp';
            end
            if vflag == 1
                vname = obj.DataObj.results.vnames(vindex,:);
            elseif vflag == 0
                vname = obj.DataObj.results.vnames(vindex+1,:);
            end;

            if nobs > 1000
                msize = 4;
            elseif nobs > 100
                msize = 8;
            elseif nobs > 50
                msize = 16;
            else
                msize = 32;
            end;
            Q0 = obj.DataObj.results.Q0;
            Q1 = obj.DataObj.results.Q1;
            Q2 = obj.DataObj.results.Q2;
            Q3 = obj.DataObj.results.Q3;
            Q4 = obj.DataObj.results.Q4;

            if obj.HFig == 0
                if obj.DataObj.results.legendmenu == 0
                    obj.HFig = figure('Position',[obj.DataObj.results.width+100 100 500 500], ... % [left bottom width height]
                        'NumberTitle','off', ...
                          'Name','Moran Scatterplot', ...
                        'Toolbar','figure');
                elseif obj.DataObj.results.legendmenu == 1
                    obj.HFig = figure('Position',[obj.DataObj.results.width+100 100 500 500], ... % [left bottom width height]
                        'NumberTitle','off', ...
                        'Name','Moran Scatterplot',...
                        'Toolbar','figure');
                end;
    
                spop = uicontrol('Style', 'popup',...
                'String', 'Select rectangle|clear selection',...
                'Position', [0 5 150 20],...
                'Callback', @(src,evnt)moranselection(obj,src,evnt));
                 chk = uicontrol('Style', 'checkbox',...
                'String', 'Map this graph',...
                'Position', [300 5 150 20],...
                'Callback', @(src,evnt)mapmoran(obj,src,evnt));
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
                obj.chx = chk;
            else
                figure(obj.HFig);
            end;

            hc = colormap(cmap);
            % hc is always 64 by 3 matrix
            incr = floor(64/nbc);
            cindex = 1:incr:64;
            cindex = cindex(1:nbc);
            hcolor = hc(cindex,:);

            map_colors = ones(nobs,3);
            for i=1:3;
                map_colors(Q0,i) = 1;
                map_colors(Q1,i) = hcolor(1,i);
                map_colors(Q2,i) = hcolor(2,i);
                map_colors(Q3,i) = hcolor(3,i);
                map_colors(Q4,i) = hcolor(4,i);
            end;
            obj.graph_color = map_colors;
            if obj.mapped ==1
                obj.DataObj.results.mapcolors = map_colors;
            end

            good = find(obj.DataObj.results.cmissing == 1);
            mycolor = zeros(nobs,3);
            mycolor(Q0,:) = ones(length(Q0),3);
            mycolor(Q1,:) = matmul(hcolor(1,:),ones(length(Q1),3));
            mycolor(Q2,:) = matmul(hcolor(2,:),ones(length(Q2),3));
            mycolor(Q3,:) = matmul(hcolor(3,:),ones(length(Q3),3));
            mycolor(Q4,:) = matmul(hcolor(4,:),ones(length(Q4),3));

            tmpvariable = cvariable1(obs_selected,vindex);
            tmpWX = WX(obs_selected,vindex);
            hscattergroup = scatter(cvariable1,tmpWX,msize,mycolor); 
            hscattergroup = scatter(cvariable1,WX(:,vindex),msize,mycolor); %get the handle of the scattergroup object created

            hpatch = get(hscattergroup,'Children');  %return the handles (array m_fig) of patches. These patchs i.e dots on the scatterplot
                                        %are scattergroup' children.

            set(hpatch(Q0),'Visible','off'); %hide all dots
            set(hpatch(Q1),'Visible','off');
            set(hpatch(Q2),'Visible','off');
            set(hpatch(Q3),'Visible','off');
            set(hpatch(Q4),'Visible','off');

            axis normal;

            %trans help to find out the right handles according to obs_selected
            trans = ones(length(obs_selected),1);
            trans = trans*(length(Q0)+length(Q1)+length(Q2)+length(Q3)+length(Q4)+1);
            trans = matsub(trans,obs_selected);

            set(hpatch(trans),'Visible','on');

            set(gcf,'UserData',hpatch);

            % try keying this only on selected observations
            % for the case of a zoom
            tmpvariable = cvariable1(obs_selected,1);
            tmpWX = WX(obs_selected,vindex);

            xbuffer = 0.5*abs(0.1*min(tmpvariable));
            ybuffer = 0.1*min(tmpWX);

            maxx = max(tmpvariable);
            minx = min(tmpvariable);
            maxy = max(tmpWX);
            miny = min(tmpWX);

            axis([1.1* minx 1.2*maxx 1.1*miny 1.2*maxy]);
            line([1.1*minx 1.2*maxx],[0 0],'Color','black');
            line([0 0],[1.1*miny 1.2*maxy],'Color','black');
            axis manual;

            if obj.DataObj.results.labels == 1
                for i=1:length(obs_selected);
                    j = obs_selected(i,1);
                    hi = text(cvariable1(j)+xbuffer,WX(j),num2str(j));
                    set(hi,'fontsize',6);
                end;
            end;
            set(obj.HFig,'Visible','on');
            obj.Patch = hpatch;
            xlabel(vname);
            mname = ['W*',vname];
            ylabel(mname);
       end
       
%   Name: updatePlot
%   Function: Update current moranplot
%   Input: obj = current Moranplot object
%   Output: Moranplot object is updated with new selected data points

       function updatePlot(obj)
           obs_selected = [];
           cmap = obj.DataObj.results.cmap;
           nbc = 4;
            vindex = obj.DataObj.results.vindex;
            vflag = obj.DataObj.results.vflag;

            cvariable1 = obj.cvariable1(:,vindex);
            obs_selected = obj.DataObj.results.obs_selected;
            nobs = length(obs_selected);
            if nobs == 0
                nobs = obj.DataObj.results.npoly; 
                temp = 1:obj.DataObj.results.npoly;
                obs_selected = temp';
            end
            if vflag == 1
                vname = obj.DataObj.results.vnames(vindex,:);
            elseif vflag == 0
                vname = obj.DataObj.results.vnames(vindex+1,:);
            end;

            if nobs > 1000
                msize = 4;
            elseif nobs > 100
                msize = 8;
            elseif nobs > 50
                msize = 16;
            else
                msize = 32;
            end;
            WX = obj.DataObj.results.WX;% we do not re-caluculate mean based on the selected map polygons
            Q0=find(cvariable1 == 0 & WX(:,vindex)== 0);
            Q1=find(cvariable1>0 & WX(:,vindex)>0);
            Q2=find(WX(:,vindex)>0 & cvariable1<= 0);
            Q3=find(WX(:,vindex)<= 0 & cvariable1<= 0);
            Q4=find(cvariable1>0 & WX(:,vindex)<= 0);

            obj.DataObj.results.Q0 = Q0;
            obj.DataObj.results.Q1 = Q1;
            obj.DataObj.results.Q2 = Q2;
            obj.DataObj.results.Q3 = Q3;
            obj.DataObj.results.Q4 = Q4;

        if obj.HFig == 0
            if obj.DataObj.results.legendmenu == 0
                obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 400 400], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Moran Scatterplot', ...
                      'MenuBar','none');
            elseif obj.DataObj.results.legendmenu == 1
                obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 400 400], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Moran Scatterplot');
            end;
    
            spop = uicontrol('Style', 'popup',...
                'String', 'Select rectangle|clear selection',...
                'Position', [0 5 150 20],...
                'Callback',  @(src,evnt)moranselection(obj,src,evnt));
            chk = uicontrol('Style', 'checkbox',...
                'String', 'Map this graph',...
                'Position', [300 5 150 20],...
                'Callback', @(src,evnt)mapmoran(obj,src,evnt));
            obj.chx = chk;

        else
            figure(obj.HFig);
        end;

        hc = colormap(cmap);
        % hc is always 64 by 3 matrix
        incr = floor(64/nbc);
        cindex = 1:incr:64;
        cindex = cindex(1:nbc);
        hcolor = hc(cindex,:);

        map_colors = ones(nobs,3);
        for i=1:3;
            map_colors(Q0,i) = 1;
            map_colors(Q1,i) = hcolor(1,i);
            map_colors(Q2,i) = hcolor(2,i);
            map_colors(Q3,i) = hcolor(3,i);
            map_colors(Q4,i) = hcolor(4,i);
        end;
        obj.graph_color = map_colors;
        if obj.mapped ==1
                obj.DataObj.results.mapcolors = map_colors;
        end

        good = find(obj.DataObj.results.cmissing == 1);
        mycolor = zeros(nobs,3);
        mycolor(Q0,:) = ones(length(Q0),3);
        mycolor(Q1,:) = matmul(hcolor(1,:),ones(length(Q1),3));
        mycolor(Q2,:) = matmul(hcolor(2,:),ones(length(Q2),3));
        mycolor(Q3,:) = matmul(hcolor(3,:),ones(length(Q3),3));
        mycolor(Q4,:) = matmul(hcolor(4,:),ones(length(Q4),3));

        % try keying this only on selected observations
        % for the case of a zoom

        tmpWX = WX(obs_selected,vindex);

        hscattergroup = scatter(cvariable1,WX(:,vindex),msize,mycolor); %get the handle of the scattergroup object created

        hpatch = get(hscattergroup,'Children');  %return the handles (array m_fig) of patches. These patchs i.e dots on the scatterplot
                                        %are scattergroup' children.

        set(hpatch(Q0),'Visible','off'); %hide all dots
        set(hpatch(Q1),'Visible','off');
        set(hpatch(Q2),'Visible','off');
        set(hpatch(Q3),'Visible','off');
        set(hpatch(Q4),'Visible','off');


        axis normal;

        %trans help to find out the right handles according to obs_selected
        trans = ones(length(obs_selected),1);
        trans = trans*(length(Q0)+length(Q1)+length(Q2)+length(Q3)+length(Q4)+1);
        trans = matsub(trans,obs_selected);

        set(hpatch(trans),'Visible','on');

        tmpvariable = cvariable1(obs_selected,1);
        xbuffer = 0.5*abs(0.1*min(tmpvariable));
        ybuffer = 0.1*min(tmpWX);

        maxx = max(tmpvariable);
        minx = min(tmpvariable);
        maxy = max(tmpWX);
        miny = min(tmpWX);

        axis([1.1* minx 1.2*maxx 1.1*miny 1.2*maxy]);

        line([1.1*minx 1.2*maxx],[0 0],'Color','black');
        line([0 0],[1.1*miny 1.2*maxy],'Color','black');
        axis manual;

        if obj.DataObj.results.labels == 1
            for i=1:length(obs_selected);
                j = obs_selected(i,1);
                hi = text(cvariable1(j)+xbuffer,WX(j,vindex),num2str(j));
                set(hi,'fontsize',6);
            end;
        end;

        set(obj.HFig,'Visible','on');
        xlabel(vname);
        mname = ['W*',vname];
        ylabel(mname);
        obj.Patch = hpatch;
  end
       
%   Name: moranselection
%   Function: Handle rectangle selection on the moran scatterplot
%   Input: obj = current Moranplot object
%           src = a source object 
%           evnt = a event object
%   Output: The information about selected data points (obs_selected) in
%   the DataSource object is updated

       function moranselection(obj,src,evnt)
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

            xlim1 = min(x);
            xlim2 = max(x);
            ylim1 = min(y);
            ylim2 = max(y);

            maxx = xlim2;
            maxy = ylim2;
            minx = xlim1;
            miny = ylim1;
            obs_selected = [];
            nobs = obj.DataObj.results.nobs;
            for i = 1:length(obj.Patch)
                j =nobs-i+1;
                x1 = get(obj.Patch(j),'XData');
                y1 = get(obj.Patch(j),'YData');
                if x1 > minx && x1 <= maxx && y1 > miny && y1 <= maxy
                    obs_selected = [obs_selected
                                               i];
                end
            end
            obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected; 
            new_variable = variable(obs_selected,:);
            new_missing = missing(obs_selected,1); 
            obj.DataObj.results.obs_selected = obs_selected;
        else 
            obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected; 
            obj.DataObj.results.obs_selected = [];
            new_variable = variable;
            new_missing = missing;
        end
        obj.DataObj.results.cvariable1 = new_variable(:,vindex); % we send down the current variable to create the legend
        obj.DataObj.results.svariable = new_variable; % contains all variables zoomed
        obj.DataObj.results.cmissing = new_missing;
        notify(obj,'graph2map');
       end 
       
       function mapmoran(obj,src,evnt)
       if (get(src,'Value') == get(src,'Max'))
        obj.mapped = 1;
       else
        obj.mapped = 0;
       end
        notify(obj,'mapping');
        end
   end
   
end 
