function poly = make_map(results)
% PURPOSE: constucts a map using structure variable returned by read_shape()
% --------------------------------------------------------------------------
% USAGE: poly = make_map(results),
% where: results is a structure variable returned by read_shape()
% --------------------------------------------------------------------------
% RETURNS: a structure variable with handles to the map polygons
%  poly(i).handles(k) = a handle to each of (nobs=npoly) polygon regions and its k-parts
%  poly(1).fig_handle = a handle to a figure containing the map
%            (use: set(poly(1).fig_handle,'Visible','on') to see the map)
% --------------------------------------------------------------------------
% NOTES:
% 1) to load and plot a map involving npoly=nobs sample data observations
% results = shape_read('myarcfile');
% poly = make_map(results);
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
poly(cnt).handles(1,1) = handles(jj);
cnt = cnt+1;
jj = jj+1;
else
 for k=1:results.nparts(cnt);
  poly(cnt).handles(1,k) = handles(jj);
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


