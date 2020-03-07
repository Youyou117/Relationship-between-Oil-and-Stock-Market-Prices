classdef am_Histogram<am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       Produce a histogram with user-specified variables and a color theme
%       for the coupled BaseMap object
% ----------------------------------------------------------------------
% see also: Representation
% ----------------------------------------------------------------------

   properties
       HV
       graph_color = [];
       chx;
   end
   properties (SetObservable)
       mapped;
   end
%   HV = a scalar indicate the direction of histogram
%           0 = horizontal, 1 = vertical

   methods
       
%   Name: Histogram
%   Function: Construction function of class 'Histogram'
%   Input: dataObj = a DataSource object to be mapped
%   Output: obj = a Histogram object 

       function obj = am_Histogram(dataObj)
           obj = obj@am_Representation(dataObj);
           obj.HFig = 0;
           obj.HV = 0;
           obj.mapped =0;
           obj.onPlot;
       end
       
%   Name: onPlot
%   Function: Produce a histogram 
%   Input: obj = current Histogram object 
%   Output: Call sub-procedure drawHistogram

       function onPlot(obj)
           obj.drawHistogram%call drawhistogram function
       end
       
%   Name: updatePlot
%   Function: Update current histogram
%   Input: obj = current Histogram object
%   Output: Call sub-procedure drawHistogram to redraw the histogram

       function updatePlot(obj)
            obj.drawHistogram%call the same drawhistogram function
       end
       
%   Name: histoselection
%   Function: Handle mouse-click on histogram bars
%   Input: obj = current Histogram object
%           src = a source object 
%           evnt = a event object
%   Output: Information about selected objects in the DataSource object is updated   

       function histoselection(obj,src,evnt)
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
            variable1 = obj.DataObj.results.variable1;
            cmissing = obj.DataObj.results.globalcmissing;
            wvariable = obj.DataObj.results.wvariable;
            vindex = obj.DataObj.results.vindex;

            val = get(src,'Value');
            obj.DataObj.results.highlightpoly1 = obj.DataObj.results.highlightpoly;
            obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected;
            if val == 1 % select bar
                k = waitforbuttonpress;
                point1 = get(gca,'CurrentPoint');  % button down detected
                pt = point1(1,1:2);    % extract x and y
                bar_h = obj.DataObj.results.bar_h; % get the handlers of bars
               
                obj.DataObj.results.highlightpoly =[];
                obj.DataObj.results.highlightclass = 0;
                for i = 1:length(bar_h)
                    x = get(bar_h(i),'XData');
                    y = get(bar_h(i),'YData');
                    if pt(1)<x(3,1) && pt(1)>x(1,1) && pt(2)>y(1,1) && pt(2)<y(3,1) %determine which bar is clicked
                        obj.DataObj.results.highlightclass = i;
                    break;
                    end             
                end

            nbc = obj.DataObj.results.nbc;
            edge=[min(variable1):(max(variable1)-min(variable1))/nbc:max(variable1)-(max(variable1)-min(variable1))/nbc];
            edge2=[edge,inf];
            bari = obj.DataObj.results.bari;
            if mod(i,2) ~=0  %if bars representing selected polygons are clicked
                i= ceil(i/2);
                i = bari(i);
            end

            for j = 1:length(obj.DataObj.results.obs_selected)
                 k = obj.DataObj.results.variable(obj.DataObj.results.obs_selected(j),vindex);%data volume is small for histogram, therefore we do not do trade off between time and space here
                if obj.DataObj.results.globalcmissing(obj.DataObj.results.obs_selected(j)) ==1
                    if k>=edge2(i) && k<edge2(i+1)   
                        obj.DataObj.results.highlightpoly = [obj.DataObj.results.highlightpoly
                                                            obj.DataObj.results.obs_selected(j)];
                    end
                end
            end
            
            else
        obj.DataObj.results.obs_selected = [];
        obj.DataObj.results.highlightpoly = [];
        obj.DataObj.results.highlightclass = 0;
        obj.DataObj.results.cmissing = cmissing;
                    end
         

        notify(obj,'graph2map'); %send out 'graph2map' events
      
           end
       
%   Name: dirselection
%   Function: Change direction of the histogram
%   Input: obj = current Histogram object
%           src = a source object 
%           evnt = a event object
%   Output: The direction of histogram is altered      

       function dirselection(obj,src,evnt)
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
            variable1 = obj.DataObj.results.variable1;
            wvariable = obj.DataObj.results.wvariable;
            vindex = obj.DataObj.results.vindex;

            val = get(src,'Value');
            if val == 1
                if(obj.HV ==0)
                    obj.HV = 1;
                end
            else
                if(obj.HV == 1)
                    obj.HV =0;
                end
            end
        obj.drawHistogram
       end
       
       %   Name: cateselection
%   Function: Handle selection of different categories
%   Input: obj = current BaseMap object
%           src = a source object 
%           evnt = a event object
%   Output: selection information in the DataSource object is updated   

       function cateselection(obj,src,evnt)
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
            obj.DataObj.results.highlightpoly1 = [];
            obj.DataObj.results.highlightpoly = [];
            obj.DataObj.results.highlightclass = 0;
            if val == 2
                nbc = 2;
            elseif val == 3
                nbc = 3;
            elseif val == 4
                nbc = 4;
            elseif val == 5
                nbc = 5;
            elseif val == 6
                nbc = 6;
            elseif val == 7
                nbc = 7;
            elseif val == 8
                nbc = 8;
            elseif val == 9
                nbc = 9;
            elseif val == 10
                nbc = 10;
            elseif val == 11
                nbc = 11;
            elseif val == 12
                nbc = 12;
            elseif val == 13
                nbc = 13;
            elseif val == 14
                nbc = 14;
            elseif val == 15
                nbc = 15;
            elseif val == 16
                nbc = 16;
            elseif val == 17
                nbc = 17;
            elseif val == 18
                nbc = 18;
            elseif val == 19
                nbc = 19;
            elseif val == 20
                nbc = 20;
            else 
                nbc = 5;
            end;

            obj.DataObj.results.nbc = nbc;
            notify(obj,'graph2map'); %send out 'graph2map' event
       end
%   Name: drawHistogram
%   Function: Produce a histogram with information in the DataSource object
%   and according to user-specified variable, category and direction
%   Input; obj = current Histogram object
%   Output: Previous Histogram is cleaned, and a new Histogram is drawn

       function drawHistogram(obj)
           cmap = obj.DataObj.results.cmap;
            nbc = obj.DataObj.results.nbc;
            
            
            missing = obj.DataObj.results.cmissing;   % missing observations in the face of a zoom

            vindex = obj.DataObj.results.vindex;
            if isempty(obj.DataObj.results.obs_selected)
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
                error('arc_mat: you need a higher screen resolution than 800x600 to use arc_mat');
            end;
            if (obj.HFig == 0)
                if obj.DataObj.results.legendmenu == 0
                        obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width+100 height+100], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Histogram', ...
                      'MenuBar','none', ...
                      'Toolbar','figure');
                elseif obj.DataObj.results.legendmenu == 1
                    obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width+100 height+100], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Histogram',...
                      'Toolbar','figure');
                end;
            else
                figure(obj.HFig);
            end;

            hc = colormap(cmap);
            clf;%clear current figure
            spop1 = uicontrol('Style', 'popup',...
                'String', 'select|clear selection',...
                'Position', [2 5 150 20],...
                'Callback',@(src,evnt)histoselection(obj,src,evnt)); %add popup menu for selection 
            spop2 = uicontrol('Style', 'popup',...
                'String', 'horizontal|vertical',...
                'Position', [160 5 100 20],...
                'Callback',@(src,evnt)dirselection(obj,src,evnt));
            cpop = uicontrol('Style', 'popup',...
                'String', '# categories|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20',...
                'Position', [300 5 100 20],...
                'Callback', @(src,evnt)cateselection(obj,src,evnt));
             chk = uicontrol('Style', 'checkbox',...
                'String', 'Map this graph',...
                'Position', [450 5 150 20],...
                'Callback', @(src,evnt)maphisto(obj,src,evnt));

             obj.DataObj.results.rect = -1.00; %initialize selection rectangle in histogram
     

            % hc is always 64 by 3 matrix
            incr = floor(64/nbc);
            cindex = 1:incr:64;
            cindex = cindex(1:nbc);
            hcolor = hc(cindex,:);

            % adjust the legend for missing values
            mindex = find(missing == 1);
            if ~isempty(mindex)
                wvariable = variable(mindex,1);
            else
                wvariable = variable;
            end;

            obj.DataObj.results.wvariable = wvariable; %local data of selected polygons

            mindex =find(obj.DataObj.results.globalcmissing ==1); %global missing represents the region with no data
            if ~isempty(mindex)
                variable1 = obj.DataObj.results.variable(mindex,vindex);
            else
                variable1 = obj.DataObj.results.variable(:,vindex);
            end
            obj.DataObj.results.variable1 = variable1; %gloabl data of all polygons

            text_h=[]; %temporary handlers
            bar_h=[];
            bari = [];

            % Trace the histogram
            edge=[min(variable1):(max(variable1)-min(variable1))/nbc:max(variable1)-(max(variable1)-min(variable1))/nbc];
            edge2=[edge,inf];

            [N,binpoints]=histc(wvariable,edge2);
            [N2,binpoints2] = histc(variable1,edge2);

            N=N(1:end-1);
            N2=N2(1:end-1);

            % form a set of nbc polygons to use with fill
            nbars = length(unique(binpoints));
            nbars2 = length(unique(binpoints2));

            edges = edge2;
            edges(end) = max(variable1);
            bari = unique(binpoints2);
            obj.DataObj.results.nbars =nbars;
            obj.DataObj.results.edges = edges;

            if ~obj.HV
                axis ij
            end
            kk=1;
            if nbars2 > 1%we need to draw rectangles and lines to represent the selected polygons
                hold on;
                k=1;
                xline(k) = edges(1); 
                yline(k) = 0;

                for i = 1:nbc
                    flag = find(bari == i, 1);
                    if ~isempty(flag)
                        j = bari(kk);
                        xc = [edges(j) edges(j) edges(j+1) edges(j+1)];
                        yc = [0 N(j) N(j) 0]; 
                        if obj.HV
                            bar_h(k) = fill(xc,yc,hcolor(j,:));%draw bars to display local data
                        else
                            bar_h(k) = fill(yc,xc,hcolor(j,:));
                        end
                        k = k+1;
                        xline(k) = edges(i);
                        yline(k) = N(i);
                        xc = [edges(j) edges(j) edges(j+1) edges(j+1)];
                        yc = [N(j) N2(j) N2(j) N(j)];
                        if obj.HV
                            bar_h(k) = fill(xc,yc,hcolor(j,:));%draw bars to display global data
                        else
                            bar_h(k) = fill(yc,xc,hcolor(j,:));
                        end
                        k = k+1;
                        xline(k) = edges(i+1);
                        yline(k) = N(i);
                        kk=kk+1;
                    else 
                        k = k+1;
                        xline(k) = edges(i);
                        yline(k) = 0;
                        k = k+1;
                        xline(k) = edges(i+1);
                        yline(k) = 0;
                    end
                end
                k = k+1;
                xline(k) = edges(i+1);
                yline(k) = 0;
    
                if size(variable1,1) ~= size(wvariable,1)
                    if obj.HV
                        plot(xline,yline,'Color',[0 1 1],'LineWidth',3);
                    else
                        plot(yline,xline,'Color',[0 1 1],'LineWidth',3);%draw line to display the upper bound of the histogram of local data
                    end
                end
                k=1;
                for i=1:nbars2;
                    j = bari(i);
                    k=k+1;

                    if  N2(j)~=0
                        if obj.HV
                            text_h(k) = text( (edges(j)+edges(j+1))/2,N2(j),num2str(N2(j)),'FontName','Tahoma','Color',[0 0 0],'FontWeight','bold');
                        else
                            text_h(k) = text( N2(j),(edges(j)+edges(j+1))/2,num2str(N2(j)),'FontName','Tahoma','Color',[0 0 0],'FontWeight','bold');
                        end  
                    else
                        text_h(k)=0;
                    end
                    k = k+1;
                end
                obj.DataObj.results.binpoints = binpoints;%update information in the datasource object
                obj.DataObj.results.binpoints2 = binpoints2;
                obj.DataObj.results.bar_h = bar_h;
                obj.DataObj.results.text_h = text_h;
                obj.DataObj.results.bari = bari;
    
                if obj.DataObj.results.highlightclass~=0%0=user have not clicked on histogram
                    pp =  obj.DataObj.results.highlightclass;
                    bartemp=obj.DataObj.results.bar_h(pp);
                    x = get(bartemp,'XData');
                    y = get(obj.DataObj.results.bar_h(pp),'YData');
                    row = size(x,1);
                    x1 = x(1:row,1);
                    y1 = y(1:row,1);
                    h = fill(x1,y1,[0.5 0.5 0.5],'EdgeColor','yellow'); %draw 'selection' rectangle on the bar clicked
                    obj.DataObj.results.rect = h; %store the handler of 'selection' rectangle        
                end
                xlabel(vname);
                hold off;

                nobs = length(obj.DataObj.results.variable);
                legend_colors = zeros(nobs,3);
                cnt = 1;

                for i=1:nobs; %determine colors for histogram bars
                    if obj.DataObj.results.globalcmissing(i) == 1
                        legend_colors(i,:) = hcolor(binpoints2(cnt,1),:);
                        cnt = cnt+1;
                    else
                        legend_colors(i,:) = [1 1 1];
                    end;
                end;
                obj.graph_color = legend_colors;
            if obj.mapped ==1
                obj.DataObj.results.mapcolors = legend_colors;
            end
            elseif nbars == 1
            % here we have only one color, so we kludge the legend and colormap
                nobs = length(obj.DataObj.results.variable);
                legend_colors = zeros(nobs,3);
                for i=1:nobs;
                    if obj.DataObj.results.globalcmissing(i) == 1
                        legend_colors(i,:) = [1 1 1];
                    else
                        legend_colors(i,:) = [1 1 1];
                    end;
                end;
                obj.graph_color = legend_colors;
            if obj.mapped ==1
                obj.DataObj.results.mapcolors = legend_colors;
            end
            figure(obj.DataObj.results.HFig);
            clf;
            Hwarning=uicontrol('Style','text','Units','Normalized','Position',[0.1,0.1,1,0.035],'Backgroundcolor',[1 1 1], ...
                'String',['Only 1 value for this variable'],'enable','inactive','FontSize',6,'HorizontalAlignment','left');
            % warning('arc_mat: you must select at least 2 polygons for the% legend');
            end
            hd = uicontextmenu; %contextmenu for sychronization choice
            obj.HEnableCm = uimenu(hd,...
                'Label','Synchronization',...
                'Checked','on',...
                'Callback',@(src,evnt)enableSyn(obj,src,evnt));
            obj.HDisableCm = uimenu(hd,...
                'Label','Non-synchronization',...
                'Checked','off',...
                'Callback',@(src,evnt)disableSyn(obj,src,evnt));
           set(obj.HFig,'uicontextmenu',hd);
           obj.chx = chk;
       end
       function maphisto(obj,src,evnt)
       if (get(src,'Value') == get(src,'Max'))
        obj.mapped = 1;
       else
        obj.mapped = 0;
       end
                notify(obj,'mapping');
        end
   end
end 
