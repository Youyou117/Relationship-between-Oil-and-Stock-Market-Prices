classdef am_LocalMoran<am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
% The LocalMoran class produces a map representing the LISA estimates -
% Local Indicator of Spatial Autocorrelation
% ----------------------------------------------------------------------
% see also: Representation
% ----------------------------------------------------------------------

   properties
       LISAMoran
       LISAMoran1
       E
       Var
   end
%   LISAMoran = a vector stores the variables under investigration
%   LISAMoran1 = a vector stores the selected subset of variables under
%   investigation
%   E = expected value in the calculation
%   Var = variance in the calculation

   methods

%   Name: LocalMoran
%   Function: Construction function of class 'LocalMoran'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a LocalMoran object 

       function obj = am_LocalMoran(dataObj)
           obj = obj@am_Representation(dataObj);
           obj.onPlot;
       end
       
%   Name: onPlot
%   Function: Call sub-procedure 'computeLocalMoran' to compute LISA
%   estimates, and call sub-procedure 'drawPlot' to draw LISA map
%   Input: obj = current LocalMoran object 
%   Output: Produce a LISA map 

       function onPlot(obj)
           obj.computeLocalMoran;
           obj.drawPlot;           
       end

%   Name: updatePlot
%   Function: Call sub-procedure 'computeLocalMoran' to compute LISA
%   estimates, and call sub-procedure 'drawPlot' to draw LISA map
%   Input: obj = current LocalMoran object 
%   Output: Produce a LISA map with user-selected subset of data            
       
       function updatePlot(obj)
            obj.computeLocalMoran;
            obj.drawPlot;
       end
       
%   Name: computeLocalMoran
%   Function: Compute LISA 
%   Input: obj = current LocalMoran object
%   Output: LISA estimates for each data points are calculated

       function computeLocalMoran(obj)
           latt = obj.DataObj.results.xc;
           long = obj.DataObj.results.yc;

           [j,W,j] = xy2cont(long,latt);
           obj.DataObj.results.W =W;
            good = find(obj.DataObj.results.globalcmissing == 1);
            bad = find(obj.DataObj.results.globalcmissing == 0);
            [nobsm,nvarsm] = size(obj.DataObj.results.variable);
            tmp = mean(obj.DataObj.results.variable(good,:)); % adjust mean for missing values
            svariable = matsub(obj.DataObj.variable,tmp); % center variables
            svariable(bad,:) = zeros(length(bad),nvarsm);
            obj.DataObj.results.svariable = svariable; % hold all transformed variables
            obj.DataObj.results.cvariable1 = svariable(:,obj.DataObj.results.vindex); %hold transformed current variable
            cvariable1 = svariable(:,obj.DataObj.results.vindex);
            WX = zeros(nobsm,nvarsm);
            WX =obj.DataObj.results.W(:,good)*svariable(good,obj.DataObj.results.vindex);


            obj.LISAMoran = [];
            m2 = 0;
            m4 = 0;
            temp =0;
            E = [];
            Var = [];
            for i = 1:nobsm
                temp = cvariable1(i)*cvariable1(i);
                m2 = m2 + temp;
                m4 = m4+temp*temp;
                E(i) = -sum(W(i,:))/(nobsm-1);
            end
            wi2 = [];
            wikh =[];

            for i = 1:nobsm
                temp = 0;
                temp1 = 0;

                for k = 1:nobsm
                    if k~=i
                        temp = temp+W(i,k)*W(i,k);
                        for h = 1:nobsm
                            if h~=i
                                temp1 = temp1 +W(i,k)*W(i,h);
                            end
                        end
                    end  
                end
                wi2 = [wi2
                     temp];
                wikh = [wikh
                    temp1];
            end
      
            m2 = m2/nobsm;
            m4 = m4/nobsm;
            b2 = m4/(m2*m2);
            obj.LISAMoran1 = [];
            for i = 1:nobsm
                obj.LISAMoran1 = [obj.LISAMoran1
                     cvariable1(i)*WX(i,1)/m2];
                Var(i) = wi2(i)*(nobsm - b2)/(nobsm-1) + wikh(i)*(2*b2-nobsm)/((nobsm-1)*(nobsm-2)) - (sum(W(i,:))*sum(W(i,:)))/((nobsm-1)*(nobsm-1));
                obj.LISAMoran = [obj.LISAMoran
                    1-norm_prb((cvariable1(i)*WX(i,1)/m2 -E(i))/(abs(sqrt(Var(i)))))]; %compute LISA statistics
            end   
       end
       
%   Name: drawPlot
%   Function: Produce a LISA map
%   Input: obj = current LocalMoran object
%   Output: 

       function drawPlot(obj)
           cmap = obj.DataObj.results.cmap;

            nbc = obj.DataObj.results.nbc; 
            variable = obj.DataObj.results.cvariable1; % we access this in case the user has 
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
                        'Name','Local Moran', ...
                        'MenuBar','none', ...
                        'Toolbar','figure');
                elseif obj.DataObj.results.legendmenu == 1
                    obj.HFig = figure('Position',[obj.DataObj.results.width+60 100 width height], ... % [left bottom width height]
                           'NumberTitle','off', ...
                        'Name','Local Moran',...
                        'Toolbar','figure');
                end;
            else
                figure(obj.HFig);
            end;
            hc = colormap(cmap);
            clf;
            obj.DataObj.results.rect = -1.00; %initialize selection rectangle in histogram
    
           mindex =find(obj.DataObj.results.globalcmissing ==1); %global missing represents the region with no data
            if length(mindex)>0
                variable1 = obj.LISAMoran(mindex,1);
            else
                variable1 = obj.LISAMoran;
            end
            obj.DataObj.results.variable1 = variable1; %gloabl data of all polygons
            wvariable = variable1(obj.DataObj.results.obs_selected,:);
            obj.DataObj.results.wvariable=wvariable;
            incr = floor(64/5);
            cindex = 1:incr:64;
            cindex = cindex(1:5);
            hcolor = hc(cindex,:);
            text_h=[];
            bar_h=[];
            bari = [];

            % Trace the histogram
            edge = [0.95,0.99,0.999,0.9999,0.99999];
            bound = edge;
            edge = 1 - (1-edge)/length(obj.LISAMoran);
            edge2=[edge,1];
            [N,binpoints]=histc(wvariable,edge2);
            [N2,binpoints2] = histc(variable1,edge2);

            % form a set of nbc polygons to use with fill
            nbars = length(unique(binpoints));
            nbars2 = length(unique(binpoints2));

            edges = edge2;
            edges(end) = max(variable1);
            bari = unique(binpoints2);
            obj.DataObj.results.nbars =nbars;
            obj.DataObj.results.edges = edges;
            edges1 = [0,1,2,3,4,5];

            kk=2;
            set(gca, 'XTick', []);   
            if nbars2 > 1
                hold on;
                k=1;
                xline(k) = edges1(1); 
                yline(k) = 0;
                for i = 1:(length(edges1)-1)
                    flag = find(bari == i, 1);
                    if isempty(flag)
                        j = bari(kk);
                        xc = [edges1(j) edges1(j) edges1(j+1) edges1(j+1)];
                        yc = [0 N(j) N(j) 0]; 
                        bar_h(k) = fill(xc,yc,hcolor(j,:)); %draw bars to display local data
                        k = k+1;
                        xline(k) = edges1(i);
                        yline(k) = N(i);
                        xc = [edges1(j) edges1(j) edges1(j+1) edges1(j+1)];
                        yc = [N(j) N2(j) N2(j) N(j)];
                        bar_h(k) = fill(xc,yc,hcolor(j,:)); %draw bars to display global data
                        pos = [(edges1(j)+edges1(j+1))/2 N2(j)];
                        text('Position',pos,'String',num2str(1-bound(j)),'Rotation',30);
                        k = k+1;
                        xline(k) = edges1(i+1);
                        yline(k) = N(i);
                        kk=kk+1;
                    else 
                        k = k+1;
                        xline(k) = edges1(i);
                        yline(k) = 0;
                        k = k+1;
                        xline(k) = edges1(i+1);
                        yline(k) = 0;
                    end
                end
                k = k+1;
                xline(k) = edges1(i+1);
                yline(k) = 0;
                if size(variable1,1) ~= size(wvariable,1)
                    plot(xline,yline,'Color',[0 1 1],'LineWidth',3); %draw line to display the upper bound of the histogram of local data
                end
                obj.DataObj.results.binpoints = binpoints;
                obj.DataObj.results.binpoints2 = binpoints2;
                obj.DataObj.results.bar_h = bar_h;
                obj.DataObj.results.text_h = text_h;
                obj.DataObj.results.bari = bari;
   
                if obj.DataObj.results.highlightclass~=0
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
            
            for i=1:nobs;
                if obj.DataObj.results.globalcmissing(i) == 1
                    if binpoints2(cnt,1)~=0
                        legend_colors(i,:) = hcolor(binpoints2(cnt,1),:);
                    else
                        legend_colors(i,:) = [1 1 1];
                    end
                    cnt = cnt+1;
                else
                    legend_colors(i,:) = [1 1 1];
                end;
            end;
            obj.DataObj.results.mapcolors = legend_colors;
      
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
            obj.DataObj.results.map_colors = legend_colors;

            %put a message in the legend figure window for the user
            figure(obj.DataObj.results.HFig);
            clf;
            Hwarning=uicontrol('Style','text','Units','Normalized','Position',[0.1,0.1,1,0.035],'Backgroundcolor',[1 1 1], ...
                 'String',['Only 1 value for this variable'],'enable','inactive','FontSize',6,'HorizontalAlignment','left');
            % warning('arc_mat: you must select at least 2 polygons for the legend');       
         end;
       end
   end
end
