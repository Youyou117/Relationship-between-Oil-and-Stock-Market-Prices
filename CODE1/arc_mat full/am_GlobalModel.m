classdef am_GlobalModel<am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
% The GlobalModel classes present the global estimates 
% of different spatial regression models, while the local estimates
% generated from selected subset of data are presented with LocalModel
% classes
% ----------------------------------------------------------------------
% see also: Representation, LocalModel
% ----------------------------------------------------------------------

   properties
       ModelName;
       graph_color;
       chx;
       mapped;
   end
%   ModelName = a string contains the name of spatial regression model

   methods%   Name: GlobalModel
%   Function: Construction function of class 'GlobalModel'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a GlobalModel object 
       


       function obj = am_GlobalModel(dataObj)
           obj = obj@am_Representation(dataObj);%call construction function of superior class 'Representation'
           obj. ModelName = ''; %initial value for modelname is null
           obj.chx = 0;
           obj.mapped = 0;
           obj.graph_color = [];
           
           y = obj.DataObj.results.variable(:,1);
           x = obj.DataObj.results.variable(:,2:end);
           
           obj.DataObj.results.yinp = y;
           obj.DataObj.results.xinp = x;
           latt = obj.DataObj.results.xc;
           long = obj.DataObj.results.yc;
           if isempty(obj.DataObj.results.W)
                [j,W,j] = xy2cont(long,latt);
                obj.DataObj.results.W =W;
           end
         
       end
       
%   Name: initialize
%   Function: virtual function to be implemented by derived model classes

       function initialize(obj)
       end
       
%   Name: callModel
%   Function: Call underlying spatial econometric functions to compute
%   spatial regression model estimates
%   Input: obj = current GlobalModel object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the model estimates

       function sresults = callModel(obj,y,x,W)
       end


%   Name: onPlot
%   Function: Produce a tabular view of model estimates
%   Input: obj = current GlobalModel object
%   Output; 
       function onPlot(obj)
           %the formats for estimation results from SAC,SAR,SEM are the
           %same, and onPlot function may be changed to display estimation
           %results from other models, for example, change the column
           %number and row number
           vnames = obj.DataObj.results.vnames;
            incr = 0.04;
            y = obj.DataObj.results.yinp;
            x = obj.DataObj.results.xinp;
            W = obj.DataObj.results.W;%spatial weight matrix
            [nobs,nvars] = size(x);
            
            sresults = obj.callModel(y,x,W);
            
            estimates = obj.getEstimate(y,x,W,sresults);
           
            obj.graph_color = obj.getColor(estimates);

            vnames1 = strvcat('Rho=','Sige=','R^2=','LogL=', 'Nobs='); 

            %roundoff digits
            for p = 1:size(sresults.tstat,1)
                sresults.tstat(p) = roundoff(sresults.tstat(p),2);
            end

            for q = 1:size(sresults.beta,1)
                sresults.beta(q) = roundoff(sresults.beta(q),4);
            end
           
           if obj.mapped ==1
                obj.DataObj.results.mapcolors = obj.graph_color;
            end
            
            
            q=q+1;
            q=q+1;
            sresults.beta(q) = roundoff(sresults.rho,4);
            q=q+1;
            sresults.beta(q) = roundoff(sresults.sige,4);
            q=q+1;
            sresults.beta(q) = roundoff(sresults.rsqr,2);
            q=q+1;
            sresults.beta(q) = roundoff(sresults.lik,4);
            q=q+1;
            sresults.beta(q) = sresults.nobs;

            rows = [size(vnames,1)-1 size(sresults.beta) size(sresults.tstat)];
            maxrow = max(rows);
            %convert estimates and varibale names into cell arrays
            data1 =cellstr(vnames);
            data2 = num2cell(sresults.beta);
            data3 = num2cell(sresults.tstat);
            data4 = cellstr(vnames1);

            %create a cell array for data to be displayed
            data = cell(maxrow,3);

            %combine cell arrays
            for i = 1:(size(data1,1)-1)
                data(i,1) = data1(i+1);
            end
            for j = 1:size(data2,1)
                data(j,2) = data2(j);
            end
            for k = 1:size(data3,1)
                data(k,3) = data3(k);
            end
            for l = 0:(size(data4,1)-1)
                data((maxrow-l),1) = data4(size(data4,1)-l);
            end

            %set properties of uitable2 
            property.colLabels = {'Variables', 'Beta','TStat'}; %column names
            property.rowNumbers = 1; %secify row number
            property.nVisibleRows = maxrow; %set all rows to be visible
            property.nVisibleCols = 3; %set all columns to be visible
            property.modal = 0; %uitable function under modal 0 will return haxes value
            property.disabledColumns = [1,2,3]; % set all columns to be unediatable
            %Properties:
            %     
            %               modal:  0 or 1 specifying modal state
            %        nVisibleRows:  int specifying number of visible rows
            %        nVisibleCols:  int specifying number of visible columns
            %           colLabels:  0, 1, or cell array of strings for column labels
            %          checkBoxes:  0, 1, or vector of int for check boxes
            %          rowNumbers:  0 or 1 for row numbers
            %
            %     disabledColumns:  [] or vector of int for preventing column editing
            %            colWidth:  1 or vector of int for indicating relative column widths
            %       highlightCell:  0 or 1  -- highlight selected cell
            %        highlightRow:  0 or 1  -- highlight row of selected cell
            %        highlightCol:  0 or 1  -- highlight column of selected cell
            %            fontsize:  int for font size
            %           precision:  s for formatting data (see fprintf)
            % 

            haxes = uitable2(data,property);
            %set figure properties

            % haxes = UITABLE('Data',data);
            obj.HFig = get(haxes,'Parent');%get figure handler from uitable's handler 
            
            set(obj.HFig,'Visible','off');
            pos = get(obj.HFig,'Position');
            pos4 = pos(4);
            pos(4) = pos(4)+25;
           
            if obj.DataObj.results.legendmenu == 0
                set(obj.HFig,'NumberTitle','off', ...
                        'Position',pos,...
                        'Name',obj.ModelName, ...
                        'MenuBar','none', ...
                        'Visible', 'off', ...
                        'BackingStore','off', ...
                        'DoubleBuffer', 'on',...
                        'Resize','on',...
                        'Toolbar','figure');
            elseif obj.DataObj.results.legendmenu == 1
                set(obj.HFig,'NumberTitle','off', ...
                        'Position',pos,...
                      'Name',obj.ModelName, ...
                      'Visible', 'off', ...
                      'BackingStore','off', ...
                      'DoubleBuffer', 'on',...
                      'Resize','on',...
                      'Toolbar','figure');
            end;
                         chk = uicontrol('Style', 'checkbox',...
                'String', 'Map this graph',...
                'Position', [0 pos4+1 150 20],...
                'Callback', @(src,evnt)mapmodel(obj,src,evnt));
            
            set(obj.HFig,'Visible','on');
            obj.chx = chk;
            
       end
       function estimates = getEstimate(obj,y,x,sresults)
       end
       function mapmodel(obj,src,evnt)
       if (get(src,'Value') == get(src,'Max'))
        obj.mapped = 1;
       else
        obj.mapped = 0;
       end
                notify(obj,'mapping');
       end
       function graph_color = getColor(obj,estimates)
            nbc = 5;
            nobs = size(estimates,1);
            hcolor = zeros(nbc,3);
            hcolor(1,:) = [0.1 0.1 0.1];
            hcolor(2,:) = [0.3 0.3 0.3];
            hcolor(3,:) = [0.5 0.5 0.5];
            hcolor(4,:) = [0.7 0.7 0.7];
            hcolor(5,:) = [0.9 0.9 0.9];
                
            edge=[min(estimates):(max(estimates)-min(estimates))/nbc:max(estimates)-(max(estimates)-min(estimates))/nbc];
            edge2=[edge,inf];

            [N,binpoints]=histc(estimates,edge2);
            [N2,binpoints2] = histc(estimates,edge2);
            
             model_colors = zeros(nobs,3);
                cnt = 1;

                for count=1:nobs; 
                    if obj.DataObj.results.globalcmissing(count) == 1
                        model_colors(count,:) = hcolor(binpoints2(cnt,1),:);
                        cnt = cnt+1;
                    else
                        model_colors(count,:) = [1 1 1];
                    end;
                end;
                graph_color = model_colors;
       end
   end
end 
