classdef am_Scatter3d<am_Representation
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
       Patch1;
       Patch2;
       spop1;
       spop2;
       spop3;
       index;
   end

%   Patch = a vector stores the handlers of scatter points

   methods
       function obj = am_Scatter3d(dataObj)
            obj = obj@am_Representation(dataObj);
           obj.initialize
           obj.HFig = 0;
           obj.index = ones(3,1);
           obj.Patch1 =[];
           obj.Patch2 = [];
           obj.onPlot;
           obj.HEnableSyn=true;
       end
       
       function initialize(obj)

       end
       
       function onPlot(obj)
            results = obj.DataObj.results;
           cmap = results.cmap;
            nbc = results.nbc;
           vflag = results.vflag; 
           nvarsm = results.nvarsm;
           vnames = results.vnames;
           variable = results.variable;
            
            
            missing = results.cmissing;   % missing observations in the face of a zoom

            vindex = results.vindex;
            
            if (size(vnames,1)<3)
                 error('arc_mat: you need at least three variables to use 3d scatterplot');
            end
           
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
                      'Name','3D Scatterplot', ...
                      'MenuBar','none', ...
                      'Toolbar','figure');
                elseif obj.DataObj.results.legendmenu == 1
                    obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width+100 height+100], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','3D Scatterplot',...
                      'Toolbar','figure');
                end;
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

     
              
            spop = uicontrol('Style', 'popup',...
                'String', 'select|clear selection',...
                'Position', [2 5 150 20],...
                'Callback',@(src,evnt)scatterselection(obj,src,evnt)); %add popup menu for selection 
            spop1 = uicontrol('Style', 'popup',...
                'String', vlist,...
                'Value',2,...
                'Position', [160 5 100 20],...
                'Callback',@(src,evnt)variableselection(obj,src,evnt));
              spop2 = uicontrol('Style', 'popup',...
                'String', vlist,...
                'Value',3,...
                'Position', [300 5 100 20],...
                'Callback',@(src,evnt)variableselection(obj,src,evnt));
              spop3 = uicontrol('Style', 'popup',...
                'String', vlist,...
                'Value',4,...
                'Position', [440 5 100 20],...
                'Callback',@(src,evnt)variableselection(obj,src,evnt));
            
            obj.index(1,1) = 1;
            obj.index(2,1) = 2;
            obj.index(3,1) = 3;
            
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
           obj.spop1 = spop1;
           obj.spop2 = spop2;
           obj.spop3 = spop3;
                
            else
                figure(obj.HFig);
            end;
            
            variable = results.variable;
            x = variable(:,obj.index(1));
            y = variable(:,obj.index(2));
            z = variable(:,obj.index(3));
            
            obs_selected1 = obj.DataObj.results.obs_selected1;
            
            handle = scatter3(x,y,z,'MarkerEdgeColor',[0.3 0.3 0.3]);          
       end
   
       function updatePlot(obj)
            results = obj.DataObj.results;
            cmap = results.cmap;
                      
            obs_selected = results.obs_selected;
             if (obj.HFig == 0)
                if obj.DataObj.results.legendmenu == 0
                        obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width+100 height+100], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','3D Scatterplot', ...
                      'MenuBar','none', ...
                      'Toolbar','figure');
                elseif obj.DataObj.results.legendmenu == 1
                    obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width+100 height+100], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','3D Scatterplot',...
                      'Toolbar','figure');
                end;
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

     
            hc = colormap(cmap);
            clf;%clear current figure
            spop = uicontrol('Style', 'popup',...
                'String', 'select|clear selection',...
                'Position', [2 5 150 20],...
                'Callback',@(src,evnt)scatterselection(obj,src,evnt)); %add popup menu for selection 
            spop1 = uicontrol('Style', 'popup',...
                'String', vlist,...
                'Value',2,...
                'Position', [160 5 100 20],...
                'Callback',@(src,evnt)variableselection(obj,src,evnt));
              spop2 = uicontrol('Style', 'popup',...
                'String', vlist,...
                'Value',3,...
                'Position', [300 5 100 20],...
                'Callback',@(src,evnt)variableselection(obj,src,evnt));
              spop3 = uicontrol('Style', 'popup',...
                'String', vlist,...
                'Value',4,...
                'Position', [440 5 100 20],...
                'Callback',@(src,evnt)variableselection(obj,src,evnt));
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
           obj.spop1 = spop1;
           obj.spop2 = spop2;
           obj.spop3 = spop3;
                
            else
                figure(obj.HFig);
            end;
            
            variable = results.variable;
            x = variable(:,obj.index(1));
            y = variable(:,obj.index(2));
            z = variable(:,obj.index(3));
                    
            obs_selected1 = obj.DataObj.results.obs_selected1;
          
           handle = scatter3(x,y,z,'MarkerEdgeColor',[0.3 0.3 0.3]);
            
            variable1 = results.variable(obs_selected,:);
            x1 = variable1(:,obj.index(1));
            y1 = variable1(:,obj.index(2));
            z1 = variable1(:,obj.index(3));
            hold on
            handle1 = scatter3(x1,y1,z1,'MarkerFaceColor',[0.6 0.6 0.6]);
            hold off        
                                             
        end


       function scatterselection(obj,src,evnt)
           
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
            temp = rotate3d;
            tst = get(temp,'Enable');%set of '3d' mode
            if strcmp(tst,'on')
                rotate3d;
            end
           val = get(src,'Value');
            variable = obj.DataObj.results.variable;
            vindex = obj.DataObj.results.vindex;
            missing = obj.DataObj.results.missing;
            obj.DataObj.results.obs_selected1 = obj.DataObj.results.obs_selected;
%             obj.DataObj.results.highlightpoly1 = [];
%             obj.DataObj.results.highlightpoly = [];
%             obj.DataObj.results.highlightclass = 0;
           if val == 1
            figure(obj.HFig);

               
               
           k = waitforbuttonpress;
           point1 = get(gca, 'CurrentPoint'); % mouse click position
           finalRect = rbbox;
           point2 = get(gca,'CurrentPoint');

ptfront = point1(1,:);
ptback = point2(2,:);

p1 = min(ptfront,ptback);             % calculate locations
offset = abs(ptfront-ptback);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
z = [p1(3) p1(3) p1(3)+offset(3) p1(3)+offset(3) p1(3)];

            xlim1 = min(x);
            xlim2 = max(x);
            ylim1 = min(y);
            ylim2 = max(y);
            zlim1 = min(z);
            zlim2 = max(z);

            maxx = xlim2;
            maxy = ylim2;
            maxz = zlim2;
            minx = xlim1;
            miny = ylim1;
            minz = zlim1;
            obs_selected = [];
            nobs = obj.DataObj.results.nobs;
            x = variable(:,obj.index(1));
            y = variable(:,obj.index(2));
            z = variable(:,obj.index(3));
            
            nobs = length(x);
            for i = 1:length(x)
                              
                x1 = x(i);
                y1 = y(i);
                z1 = z(i);
                
                if x1 > minx && x1 <= maxx && y1 > miny && y1 <= maxy && z1>minz && z1<=maxz
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
       function variableselection(obj,src,evnt)
           obj.index(1) = get(obj.spop1,'Value') - 1;
           obj.index(2) = get(obj.spop2,'Value') - 1;
           obj.index(3) = get(obj.spop3,'Value') - 1;
            notify(obj,'graph2map');
       end
       
      
   end
   
end 
