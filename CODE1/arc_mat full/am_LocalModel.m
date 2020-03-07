classdef am_LocalModel<am_Representation
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
% The LocalModel classes present the local estimates 
% generated from selected subset of data while the global estimates of
% different spatial regression models are presented with GlobalModel
% classes
% ----------------------------------------------------------------------
% see also: Representation, GlobalModel
% ----------------------------------------------------------------------
   properties
        ModelName
        graph_color;
        chx;
        mapped;
   end
%   ModelName = a string contains the name of spatial regression model
   methods

%   Name: LocalModel
%   Function: Construction function of class 'LocalModel'
%   Input: dataObj = a DataSource object to be modeled
%   Output: obj = a LocalModel object 

       function obj = am_LocalModel(dataObj)
           obj = obj@am_Representation(dataObj); %call construction function of superior class ' representation'
           obj.ModelName = '';
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
           end%compute spatial weight matrix from lattitude and longtitude, users may define their own W matrix to relace this line
           %for example, users may use 'make_neighborsw(latt,long,n)'
           %function to create n-neighbors W matrix; 
       end

%   Name: callModel
%   Function: Call underlying spatial econometric functions to compute
%   spatial regression model estimates
%   Input: obj = current LocalModel object
%           y = a vector stores dependent variable
%           x = a vector stores explanatroy variables
%           W = a matrix stores spatial contiguity matrix 
%   Output: sresults = a strucutre variable stores the model estimates
        
       function sresults = callModel(obj,y,x,W)
       end

%   Name: initialize
%   Function: virtual function to be implemented by derived model classes

       function initialize(obj)
       end

%   Name: updatePlot
%   Function: Since local models are always computed after user selection of
%   map polygons, we do not present LocalModel results with 'onPlot' function, and use 
%   'updatePlot' to present LocalModel results when selections are made. 
%   Input: obj = current LocalModel object
%   Output; Produce a tabular view of local model estimates

       function updatePlot(obj)
           if obj.HFig ~=0
               set(obj.HFig,'Visible','off');
               delete(obj.HFig);
               obj.HFig = 0;
           end
           vnames = obj.DataObj.results.vnames;
            incr = 0.04;

            
            % do local sub-model
            obs_selected = obj.DataObj.results.obs_selected;
            nobs = length(obs_selected);
            if nobs == 0
                nobs = obj.DataObj.results.npoly; 
                temp = 1:obj.DataObj.results.npoly;
                obs_selected = temp';
                return;
            end
            y = obj.DataObj.results.yinp(obs_selected,1);
            x = obj.DataObj.results.xinp(obs_selected,:);
            W = obj.DataObj.results.W(obs_selected,obs_selected);
            [nobs,nvars] = size(x);
            sresults = obj.callModel(y,x,W);
            
            estimates = obj.getEstimate(y,x,W,sresults);
           
            obj.graph_color = obj.getColor(estimates);
            
            vnames = strvcat(vnames,' ','Rho=','Sige=','R^2=','LogL=', 'Nobs=');
            % =============== local estimation results presented here ===============================================================================================
            %compute the row number in UITable
            rows = [size(vnames,1)-1 size(sresults.beta) size(sresults.tstat)];
            maxrow = max(rows);

            %roundoff digits
            for p = 1:size(sresults.tstat,1)
                sresults.tstat(p) = roundoff(sresults.tstat(p),2);
            end
            for q = 1:size(sresults.beta,1)
                sresults.beta(q) = roundoff(sresults.beta(q),4);
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

            %convert estimates and varibale names into cell arrays
            data1 =cellstr(vnames);
            data2 = num2cell(sresults.beta);
            data3 = num2cell(sresults.tstat);

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

            %set properties of uitable2 
            property.colLabels = {'Variables', 'Beta','TStat'};
            property.rowNumbers = 1;
            property.nVisibleRows = maxrow;
            property.nVisibleCols = 3;
            property.modal = 0;
            property.disabledColumns = [1,2,3];
            %Properties of uitable:
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

            %get figure handler
            haxes = uitable2(data,property);
            obj.HFig = get(haxes,'Parent');

            %shift figure for local estimates a littlt bit 
            position = get(obj.HFig,'Position');
            position(1) = position(1) - 30;
            position(1) = position(2) - 30;
            pos4 = position(4);
            position(4) = position(4)+25;

            %set figure properties
            if obj.DataObj.results.legendmenu == 0
                set(obj.HFig,'Position', position, ...
                      'NumberTitle','off', ...
                      'Name',obj.ModelName, ...
                      'MenuBar','none', ...
                      'Visible', 'off', ...
                      'BackingStore','off', ...
                      'DoubleBuffer', 'on', ...
                      'Resize','on',...
                      'Toolbar','figure');
            elseif obj.DataObj.results.legendmenu == 1
                set(obj.HFig,'Position', position, ...
                      'NumberTitle','off', ...
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
            obs_selected = obj.DataObj.results.obs_selected;
            nobtemp = length(obs_selected);
            if nobtemp == 0
                nobtemp = obj.DataObj.results.npoly; 
                temp = 1:obj.DataObj.results.npoly;
                obs_selected = temp';
            end
            nobs = size(obj.DataObj.results.yinp,1);
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
            
            model_colors = ones(nobs,3);
            nobs1 = size(estimates);
                for count=1:nobs1; 
                      model_colors(obs_selected(count),:) = hcolor(binpoints2(count,1),:);
                end;
                 
                graph_color = model_colors;
       end
   end
end 
