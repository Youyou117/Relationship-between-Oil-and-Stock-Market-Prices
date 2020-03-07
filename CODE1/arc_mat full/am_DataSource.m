classdef am_DataSource<handle
% ----------------------------------------------------------------------
% Created by;
%       Xingjian Liu & James. P. LeSage, 2009
%       Texas State University-San Marcos
%       spatial-econometrics.com   
%----------------------------------------------------------------------
% Usage:
%       DataSource object stores spatial and non-spatial information in 
%       the analysis.
% ----------------------------------------------------------------------
% see also: 
% ----------------------------------------------------------------------

   properties(SetObservable = true)
       results
       variable
       option
   end
%   variable = a variable vector (nobs x 1) or matrix (nobs x nvars)
%   results = a structure variable returned by shape_read()
%   options  = a structure variable with function options
%   These three variables are the same as those three variables used to
%   construct mapping GUIs in the previous version of Arc_Mat

   methods
%   Name: DataSource
%   Function: Construction function of class 'DataSource'
%   Input: variable = a variable vector (nobs x 1) or matrix (nobs x nvars)
%           results = a structure variable returned by shape_read()
%           options  = a structure variable with function options
%   Output: obj = a DataSource object
       function obj = am_DataSource(variable, results,options) 
           obj.results = results; % Assign property values
           obj.variable = variable;
           obj.option = options;
           [nobsm,nvarsm] = size(variable);
           missing_vec = ones(nobsm,1);
           % some error checking
           if nargin == 3 % user-defined options input
               if ~isstruct(options)
                   error('arc_histmap: must supply options as a structure variable');
               end;
           fields = fieldnames(options);
           nf = length(fields);
           nbc = 5;  % defaults
           cmap = 'hsv';
           vnames =  'Variable';
           for i=1:nvarsm;
                vnames = strvcat(vnames,['variable',num2str(i)]);
           end;
           vflag = 0;
           mapmenu = 0; legendmenu = 0;
           for i=1:nf
                if strcmp(fields{i},'nbc')
                    nbc = options.nbc; 
                elseif strcmp(fields{i},'cmap')
                    cmap = options.cmap;
                elseif strcmp(fields{i},'vnames')
                     vnames = options.vnames;
%                     [nchk,junk] = size(vnames);
%                     if nchk ~= nvarsm
%                         error(' wrong number of variable names');
%                     end;
                    vflag = 1;
                elseif strcmp(fields{i},'mapmenu')
                    mapmenu = options.mapmenu;  
                elseif strcmp(fields{i},'legendmenu');
                    legendmenu = options.legendmenu;
                elseif strcmp(fields{i},'missing');
                    missing_vec = options.missing;
                end;
            end;
            elseif nargin == 2 % set default options
                nbc = 5;
                cmap = 'hsv';
                vnames =  'Variable';
                for i=1:nvarsm;
                    vnames = strvcat(vnames,['variable',num2str(i)]);
                end;
                vflag = 0;
                mapmenu = 0; 
                legendmenu = 0;
           else
            error('arc_mat: Wrong # of input arguments');
        end;

        obj.results.labels = 1;%labels on plots
        obj.results.nbc = nbc; %number of categories on the histogram
        obj.results.nbcm = 4;
        obj.results.cmap = cmap;
        obj.results.vnames = vnames;%variable names
        obj.results.mapmenu = mapmenu;%menu toolbar
        obj.results.legendmenu = legendmenu; %plot legend
        obj.results.vindex = 1;%current selected variable index
        obj.results.legend_fig = 0;
    	obj.results.mapfigure = 0;
        obj.results.mapfigure1 = 0;
        obj.results.variable = variable;
        [nobsm,nvarsm] = size(variable);
        obj.results.nvarsm = nvarsm;%number of vaiables
        obj.results.vflag = vflag;
        mindex = find(missing_vec == 1);
        obj.results.svariable = variable; % holds all variables zoomed on the map
        obj.results.cvariable = variable(:,obj.results.vindex); % holds current variable selection
        obj.results.missing = missing_vec;%missing values
        obj.results.cmissing = missing_vec;
        obj.results.globalmissing = missing_vec;%hold the initial missing values
        obj.results.globalcmissing = missing_vec;
        obj.results.obs_selected = [];%currently selection
        obj.results.obs_selected1 = [];%previsously selection
        obj.results.hpatch = -1;%reserved variable for haptching
        obj.results.rect =-1;
        obj.results.mapcolors =[];%map colors
        obj.results.Hwarning = 0;
        obj.results.W = [];
        obj.results.bar_h = 0;% bars on the histogram
        obj.results.text_h = 0;% texts on the histogram
        obj.results.highlightclass = 0;%highlighted category on the histogram
        obj.results.highlightpoly1 = [];%currently selected bar on the histogram
        obj.results.highlightpoly = [];%previously selected bar on the histogram

    
        svec = get(0,'ScreenSize');%determine screen size
        if svec(3) > 1300
            width = 800; height = 800;
        elseif svec(3) > 1000
            width = 650; height = 650;
        elseif svec(3) == 800
            error('arc_mat: you need a higher screen resolution to use arc_map');
        end;
        obj.results.width = width;
        obj.results.height = height;
     end 
   end
end 
