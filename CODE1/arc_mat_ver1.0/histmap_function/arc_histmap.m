function arc_histmap(variable,results,options)
% PURPOSE: produce a map with histogram legend using ArcView shape files
% ----------------------------------------------------------------------
% USAGE: arc_histmap(variable,results,options)
% where: variable = a variable vector (nobs x 1) or matrix (nobs x nvars)
%        results  = a structure variable returned by shape_read()
%        options  = a structure variable with function options
%        options.vnames     = a string with the variable name or names
%        e.g. vnames = strvcat('constant','pinstruction','pbuilding','padminist');
%        options.cmap       = a colormap string, e.g., 'hsv','jet' (default 'hsv')
%                             (see colormap menu in the gui for various options)
%        options.nbc        = # of categories (default = 5)
%                             (see #categories menu in the gui for various options)
%        options.missing    = a vector of 0's and 1's for missing and non-missing observations
%                             0 = missing, 1 = non-missing 
%                             (produces white map polygons for missing values)
%        options.mapmenu    = 1, for a menubar, 0 = default, no menubar
%                          (useful for printing, editing the map)
%        options.legendmenu = 1, for a menubar, 0 = default, no menubar
%                          (useful for printing, editing the map)
% ----------------------------------------------------------------------
% RETURNS: a graphical user interface with the map and histogram legend
% as well as menus for: colormap,zoom map,#categories,quit
% NOTE: a right-mouse-click on the map polygons presents the data value
%       associated with the map region on which you clicked
% ----------------------------------------------------------------------
% see also: shape_read(), make_map()
% ----------------------------------------------------------------------

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
        [nchk,junk] = size(vnames);
        if nchk ~= nvarsm
        error('arc_histmap: wrong number of variable names');
        end;
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
    mapmenu = 0; legendmenu = 0;
else
error('arc_histmap: Wrong # of input arguments');
end;

results.nbc = nbc;
results.cmap = cmap;
results.vnames = vnames;
results.mapmenu = mapmenu;
results.legendmenu = legendmenu;
results.vindex = 1;
results.legend_fig = 0;
results.variable = variable;
[nobsm,nvarsm] = size(variable);
results.nvarsm = nvarsm;
results.vflag = vflag;
mindex = find(missing_vec == 1);
results.svariable = variable; % holds all variables zoomed on the map
results.cvariable = variable(:,results.vindex); % holds current variable selection
results.missing = missing_vec;
results.cmissing = missing_vec;

% figure out size of map windows
svec = get(0,'ScreenSize');
if svec(3) > 1300
width = 800; height = 800;
elseif svec(3) > 1000
width = 650; height = 650;
elseif svec(3) == 800
error('arc_histmap: you need a higher screen resolution to use arc_map');
end;
results.width = width;
results.height = height;

% construct a legend for the map
% uses results.cvariable
results = am_histlegend(results);

map_colors = results.map_colors;

poly = make_map(results,variable);

if vflag == 1
mname = ['Map of ',vnames(results.vindex,:)];
elseif vflag == 0
mname = ['Map of ',vnames(2,:)];
end;

if mapmenu == 0
set(poly(1).fig_handle, ...
            'Position',[50 100 width height], ... % [left bottom width height]
            'MenuBar','none', ...
            'NumberTitle','off', ...
            'Name',mname);
elseif mapmenu == 1
set(poly(1).fig_handle, ...
            'Position',[50 100 width height], ...
            'NumberTitle','off', ...
            'Name',mname);
end;


axis equal;
thandles = zeros(results.npoly,1);

hold on;
for i=1:results.npoly;
 for k=1:results.nparts(i);
            set(poly(i).handles(k),'FaceColor',map_colors(i,:),'Visible','on');
 end;
end;
hold off;


hpop = uicontrol('Style', 'popup',...
       'String', 'Colormap|hsv|hot|cool|gray|bone|copper|jet|pink|white|autumn|spring|winter|summer',...
       'Position', [0 5 75 20],...
       'Callback', 'am_histcmap');
%       'TooltipString','Select a color scheme for the map');
% [left bottom width height]
spop = uicontrol('Style', 'popup',...
       'String', 'Zoom map|select rectangle|zoom out',...
       'Position', [76 5 100 20],...
       'Callback', 'am_histmap');
%     'TooltipString','Zoom in on a region of the map');
cpop = uicontrol('Style', 'popup',...
       'String', '# categories|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20',...
       'Position', [176 5 125 20],...
       'Callback', 'am_histncat');
%       'TooltipString','Set the # of histogram categories');

vlist = 'Variable|';
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
       'Position', [300 5 150 20],...
       'Callback', 'am_histvariable');
%       'TooltipString','Select a variable to map');
   
qpop = uicontrol('Style', 'popup',...
       'String', 'Quit|Close/Exit',...
       'Position', [551 5 75 20],...
       'Callback', 'am_histquit');
%       'TooltipString','Close the map and histogram windows');

% place a bunch of stuff into the results structure
% for use by the sub-functions
results.hpop = hpop;
results.spop = spop;
results.cpop = cpop;
results.vpop = vpop;
results.qpop = qpop;

set(poly(1).fig_handle,'Visible','on');
set(poly(1).fig_handle,'UserData',poly);

guidata(hpop, results);
guidata(cpop, results);
guidata(spop, results);
guidata(vpop, results);
guidata(qpop, results);

% bring legend figure back to the front
figure(results.legend_fig);

% --------------------------------------------------------------------
% support functions 
% --------------------------------------------------------------------
% am_histmap, am_histcmap, am_histncat, am_histvariable, am_getinfo, am_histquit
% --------------------------------------------------------------------


function poly = make_map(results,data)
% PURPOSE: constucts a map using structure variable returned by read_shape()
% --------------------------------------------------------------------------
% USAGE: poly = make_map(results,data),
% where: results is a structure variable returned by read_shape()
%        data is the data matrix to be mapped (an input argument to arc_map
%        functions)
% --------------------------------------------------------------------------
% RETURNS: a structure variable with handles to the map polygons
%  poly(i).handles(k) = a handle to each of (nobs=npoly) polygon regions and its k-parts
%                       (these are patch objects)
%  poly(1).fig_handle = a handle to a figure containing the map
%            (use: set(poly(1).fig_handle,'Visible','on') to see the map)
% Also sets the UserData field for the patch objects to contain the
% associated row of data for this polygon
%
% --------------------------------------------------------------------------
% NOTES:
% 1) to load and plot a map involving npoly=nobs sample data observations
% results = shape_read('myarcfile');
% poly = make_map(results,results.data);
% set(poly(1).fig_handle,'Visible','on');
% 2) Also you can do the following:
% figure(1);
% hold on;
% h = []; % handles to each polygon
% for i=1:results.npoly;
%  for k=1:results.nparts(i);
%	mypoly.faces = get(poly(i).handles(k),'Faces');
%	mypoly.vertices = get(poly(i).handles(k),'Vertices');
%	h(i,k) = patch(mypoly);
% end;
% end;
% 3) to set the facecolor of the polygons, using the handles h
% for i=1:npoly;
%  for k=1:results(i).nparts;
%  set(h(i,k),'FaceColor',[0 1 1]);
%  end;
% end;

x = results.x;
y = results.y;

poly(1).fig_handle = figure('Visible','off');
handles = polyplot(x,y,'fill',[0 0 0]);

 % Process chunks separated by NaN .................
in = [0; find(isnan(x)); length(x)+1];
n = length(in)-1;
cnt = 1;
jj = 1;
while (jj <= n)
  ii = in(jj)+1:in(jj+1)-1;
  ii = [ii ii(1)];
  xx = x(ii); yy = y(ii);
if results.nparts(cnt) == 1
poly(cnt).handles(1) = handles(jj);
cmenu = uicontextmenu;
set(poly(cnt).handles(1),'UIContextMenu',cmenu);
    hi = uimenu(cmenu,'Label',[num2str(cnt) ') ' num2str(data(cnt,1))]); 
    set(poly(cnt).handles(1),'UserData',hi);
cnt = cnt+1;
jj = jj+1;

else    
 for k=1:results.nparts(cnt);
  poly(cnt).handles(k) = handles(jj);
  cmenu = uicontextmenu;
    set(poly(cnt).handles(k),'UIContextMenu',cmenu);
    hi = uimenu(cmenu,'Label',[num2str(cnt) ') ' num2str(data(cnt,1))]); 
    set(poly(cnt).handles(k),'UserData',hi);
  jj = jj+1;
 end;
cnt = cnt+1;
end;
end;


function [handles] = polyplot(x,y,a1,a2)
% POLYPLOT Plotting or filling polygons.
%	L = POLYPLOT(X,Y) plots polygon(s)
%	concatenated into coordinate vectors X, Y.
%	If X, Y consist of coordinates of several
%	polygons they must be separated by NaN:
%	X = [X1 NaN X2 NaN X3 ...]. In this case each
%	polygon is "closed" and plotted separately.
%	L is a vector of handles of lines defining
%	polygon boundaries, one handle per line.
%	L = POLYPLOT(X,Y,C) also specifies line color.
%	C can be a letter such as 'w', 'y', 'c', etc.,
%	a 1 by 3 vector in RGB format or a string of 
%	such letters, like 'rgby' or n by 3 matrix.
%	In the latter case this string or matrix plays the
%	role of color order for succession of polygons.
%	If the number of polygons is more than number of
%	colors, colors are cyclically repeated.
%
%	P = POLYPLOT(X,Y,'fill',C) fills polygons
%	creating a patch rather than a line and returns
%	patch handles P.
%
%	This program can also fill non-simply connected
%	polygons, such as ones with holes. It checks
%	the direction of each polygons separated by
%	NaN in coordinate vectors. If the contour is
%	clock-wise (signed area is negative) then it
%	interprets such a polygon as a "hole" and fills
%	it with the background color.

%  Copyright (c) 1995 by Kirill K. Pankratov,
%       kirill@plume.mit.edu.
%       06/25/95, 09/07/95  

%  May call AREA function.

 % Handle input ....................................
is_patch = 0;
clr = get(gca,'colororder');
if nargin>2
  lm = min(length(a1),4);
  names = ['fill '; 'patch'];
  is_patch = all(a1(1:lm)==names(1,1:lm));
  is_patch = is_patch | all(a1(1:lm)==names(2,1:lm));

  if is_patch
    if nargin>3, clr = a2; end
  else
    clr = a1;
  end
end
if isstr(clr), clr=clr(:); end
nclr = size(clr,1);
x = x(:); y = y(:);

 % Check hold state ............
if ~ishold, newplot, end

 % Setup a call ................
if is_patch
  call = 'patch';
  cpn = 'facecolor';
else
  call = 'line';
  cpn = 'color';
end
% call = ['p(jj)=' call '(''xdata'',xx,''ydata'',yy);'];
% call
 % Get color for "holes" polygons ..................
clrh = get(gca,'color');
if strcmp(clrh,'none'), clrh = get(gcf,'color'); end 

 % Process chunks separated by NaN .................
in = [0; find(isnan(x)); length(x)+1];
n = length(in)-1;
for jj=1:n
  ii = in(jj)+1:in(jj+1)-1;
  ii = [ii ii(1)];
  xx = x(ii); yy = y(ii);

  % Check area
  a(jj) = area(xx,yy);

  % Create the patch or line
  handles(jj) = patch(xx,yy,[1 0 0]);
  ic = rem(jj-1,nclr)+1;
  set(handles(jj),cpn,clr(ic,:))
end

 % If non-simply-connected patch, fill "holes" with 
 % background color ...............................
if is_patch & n>1
  % Find which polygons are inside which
  holes = find(a<0);
  % Set color
  set(handles(holes),'FaceColor',clrh)

end


function  a = area(x,y)
% AREA  Area of a planar polygon.
%	AREA(X,Y) Calculates the area of a 2-dimensional
%	polygon formed by vertices with coordinate vectors
%	X and Y. The result is direction-sensitive: the
%	area is positive if the bounding contour is counter-
%	clockwise and negative if it is clockwise.
%
%	See also TRAPZ.

%  Copyright (c) 1995 by Kirill K. Pankratov,
%	kirill@plume.mit.edu.
%	04/20/94, 05/20/95  

 % Make polygon closed .............
x = [x(:); x(1)];
y = [y(:); y(1)];

 % Calculate contour integral Int -y*dx  (same as Int x*dy).
lx = length(x);
a = -(x(2:lx)-x(1:lx-1))'*(y(1:lx-1)+y(2:lx))/2;


