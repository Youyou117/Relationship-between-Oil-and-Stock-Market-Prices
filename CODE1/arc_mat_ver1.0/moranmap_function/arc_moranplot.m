function arc_moranplot(variable,W,results,options)
% PURPOSE: produce a map with moran scatterplot using ArcView shape files
% ----------------------------------------------------------------------
% USAGE: arc_moranplot(variable,W,results,options)
% where: variable = a variable vector (nobs x 1) or matrix (nobs x nvars)
%               W = a spatial weight matrix
%         results = a structure variable returned by shape_read()
%        options  = a structure variable with function options
%        options.vnames     = a string with the variable name or names
%        e.g. vnames = strvcat('constant','pinstruction','pbuilding','padminist');
%        options.labels     = 1 for labels, 0 = default, no labels
%        options.missing    = a vector of 0's and 1's for missing and non-missing observations
%                             0 = missing, 1 = non-missing 
%                             (produces white map polygons for missing values)
%        options.mapmenu    = 1, for a menubar, 0 = default, no menubar
%                          (useful for printing, editing the map)
%        options.legendmenu = 1, for a menubar, 0 = default, no menubar
%                          (useful for printing, editing the map)
% ----------------------------------------------------------------------
% RETURNS: a graphical user interface with the map and moran scatterplot
% as well as menus for: zoom map,quit
% ----------------------------------------------------------------------
% see also: shape_read()
% ----------------------------------------------------------------------

% some error checking
if nargin == 4 % user-defined options input
	if ~isstruct(options)
        error('arc_moranplot: must supply options as a structure variable');
	end;
[nobsm,nvarsm] = size(variable);
missing_vec = ones(nobsm,1);

	fields = fieldnames(options);
	nf = length(fields);
	nbc = 4;  % defaults
	cmap = 'hsv';
	vnames =  'Variable';
    for i=1:nvarsm;
    vnames = strvcat(vnames,['variable',num2str(i)]);
    end;
    vflag = 0;
    labels = 0;
    mapmenu = 0; legendmenu = 0;

 for i=1:nf
    if strcmp(fields{i},'vnames');
        vnames = options.vnames;
        [nchk,junk] = size(vnames);
        if nchk ~= nvarsm;
        error('arc_histmap: wrong number of variable names');
        end;
        vflag = 1;
    elseif strcmp(fields{i},'labels');
        labels = options.labels;
    elseif strcmp(fields{i},'mapmenu')
        mapmenu = options.mapmenu;  
    elseif strcmp(fields{i},'legendmenu');
        legendmenu = options.legendmenu;
    elseif strcmp(fields{i},'missing');
        missing_vec = options.missing;
    end;
 end;

elseif nargin == 3 % set default options
[nobsm,nvarsm] = size(variable);
	vnames =  'Variable';
    for i=1:nvarsm;
    vnames = strvcat(vnames,['variable',num2str(i)]);
    end;
    vflag = 0;

	nbc = 4;
	cmap = 'hsv';
    mapmenu = 0; legendmenu = 0;
    labels = 0;
else
error('arc_moranplot: Wrong # of input arguments');
end;

% some error checking
ind = isnan(variable);
if sum(ind) > 0
error('arc_moranplot: variable vector contains NaN values');
end;
ind = isinf(variable);
if sum(ind) > 0
error('arc_moranplot: variable vector contains Inf values');
end;

[n1,n2] = size(W);
if n1 ~= results.npoly;
    error('arc_moranplot: wrong size W-matrix');
elseif n1~= n2;
    error('arc_moranplot: W-matrix is not square');
end;

[n1,n2] = size(variable);
if n1~= results.npoly;
    error('arc_moranplot: wrong size data matrix on input');
elseif n2~= nvarsm;
    error('arc_moranplot: wrong # of variable names on input');
end;

[n1,n2] = size(missing_vec);
if n1~= results.npoly;
    error('arc_moranplot: wrong size missing values vector on input');
elseif n2~= 1;
    error('arc_moranplot: wrong size missing values vector on input');
end;
    
  

results.nbc = nbc;
results.cmap = cmap;
results.moran_fig = 0;
results.variable = variable;
results.mapmenu = mapmenu;
results.legendmenu = legendmenu;
results.vnames = vnames;
results.labels = labels;
vindex = 1;
results.vindex = vindex;
results.vflag = vflag;
results.missing = missing_vec;
results.cmissing = missing_vec;

tmp = 1:results.npoly;
results.obs_selected = tmp';

% we compute the moranplot information only once
% if the user zooms in, we simply report a subset of this information

good = find(results.cmissing == 1);
bad = find(results.cmissing == 0);

tmp = mean(variable(good,:)); % adjust mean for missing values
svariable = matsub(variable,tmp); % center variables
svariable(bad,:) = zeros(length(bad),nvarsm);
results.svariable = svariable; % hold all variables
results.cvariable = svariable(:,results.vindex); % holds current variable selection
cvariable = results.cvariable;
WX = zeros(nobsm,nvarsm);
WX = W(:,good)*svariable(good,:);

% define the quadrants
Q0=find(cvariable == 0 & WX(:,vindex)== 0);
Q1=find(cvariable>0 & WX(:,vindex)>0);
Q2=find(WX(:,vindex)>0 & cvariable<= 0);
Q3=find(WX(:,vindex)<= 0 & cvariable<= 0);
Q4=find(cvariable>0 & WX(:,vindex)<= 0);
%%%%%%%%%%%%%%%%%%%%%%%

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
results.heigh = height;

% construct a legend for the map
results.Q0 = Q0;
results.Q1 = Q1;
results.Q2 = Q2;
results.Q3 = Q3;
results.Q4 = Q4;
results.WX = WX;
results = am_moran(results);

map_colors = results.map_colors;

poly = make_map(results,svariable);

if vflag == 1
mname = ['Map of ',vnames(results.vindex,:)];
elseif vflag == 0
mname = ['Map of ',vnames(vindex+1,:)];
end;

if mapmenu == 0
set(poly(1).fig_handle, ...
            'Position',[50 50 width height], ...
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


  if results.labels == 1	
	cnt1 = 1;
	cnt2 = 1;
	cnt3 = 1;
	cnt4 = 1;
  end;

hold on;
thandles = zeros(results.npoly,1);
for i=1:results.npoly;
 for k=1:results.nparts(i);
    set(poly(i).handles(k),'FaceColor',map_colors(i,:),'Visible','on');
 end;
  if results.labels == 1
	thandles(i,1) = text(results.xc(i),results.yc(i),num2str(i));
	set(thandles(i,1),'fontsize',6);
  end;
end;
hold off;

% [left bottom width height]
spop = uicontrol('Style', 'popup',...
       'String', 'Zoom map|select rectangle|zoom out',...
       'Position', [0 5 150 20],...
       'Callback', 'am_selectmap');

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
       'Position', [151 5 150 20],...
       'Callback', 'am_moranvariable');

qpop = uicontrol('Style', 'popup',...
       'String', 'Quit|Close/Exit',...
       'Position', [501 5 100 20],...
       'Callback', 'am_moranquit');

% place a bunch of stuff into the results structure
% for use by the sub-functions
results.spop = spop;
results.qpop = qpop;
results.vpop = vpop;

results.thandles = thandles;


set(poly(1).fig_handle,'Visible','on');
set(poly(1).fig_handle,'UserData',poly);

guidata(spop, results);
guidata(qpop, results);
guidata(vpop, results);

% --------------------------------------------------------------------
% support functions 
% --------------------------------------------------------------------
% am_moran, am_selectmap, am_moranvariable, am_quit
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
end; % end of if/else
end; % end of while


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


