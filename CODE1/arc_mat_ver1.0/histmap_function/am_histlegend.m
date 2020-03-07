function results = am_histlegend(results)
% PURPOSE: compute legend figure for histmap
% returns an n x 3 matrix of colors that can be used to label the map polygons
% and a handle to the legend figure

cmap = results.cmap;
nbc = results.nbc;
variable = results.cvariable; % we access this in case the user has 
                              % changed to a subset of the sample data 
missing = results.cmissing;   % missing observations in the face of a zoom

vindex = results.vindex;
vflag = results.vflag;
vnames = results.vnames;
if vflag == 1
vname = vnames(vindex,:);
elseif vflag == 0
vname = vnames(vindex+1,:);
end;

svec = get(0,'ScreenSize');
if svec(3) > 1300
width = 400; height = 400;
elseif svec(3) > 1000
width = 400; height = 400;
elseif svec(3) == 800
error('arc_histmap: you need a higher screen resolution than 800x600 to use arc_map');
end;
if (results.legend_fig == 0)
	if results.legendmenu == 0
	legend_fig = figure('Position',[results.width+60 100 width height], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Map Legend', ...
                      'MenuBar','none');
	elseif results.legendmenu == 1
	legend_fig = figure('Position',[results.width+60 100 width height], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Map Legend');
	end;
	results.legend_fig = legend_fig;
else
figure(results.legend_fig);
end;
hc = colormap(cmap);

% hc is always 64 by 3 matrix
incr = floor(64/nbc);
cindex = 1:incr:64;
cindex = cindex(1:nbc);
hcolor = hc(cindex,:);

% adjust the legend for missing values
mindex = find(missing == 1);
if length(mindex) > 0
wvariable = variable(mindex,1);
else
wvariable = variable;
end;

% Trace the histogram
edge=[min(wvariable):(max(wvariable)-min(wvariable))/nbc:max(wvariable)-(max(wvariable)-min(wvariable))/nbc];
edge2=[edge,inf];
[N,binpoints]=histc(wvariable,edge2);
N=N(1:end-1);
% form a set of nbc polygons to use with fill
nbars = length(unique(binpoints));
edges = edge2;
edges(end) = max(wvariable);
bari = unique(binpoints);
if nbars > 1
clf;
hold on;
	for i=1:nbars;
	j = bari(i);
	xc = [edges(j) edges(j) edges(j+1) edges(j+1)];
	yc = [0 N(j) N(j) 0];
	bar_h(i) = fill(xc,yc,hcolor(j,:));
	end;
hold off;
xlabel(vname);
nobs = length(variable);
legend_colors = zeros(nobs,3);
cnt = 1;
for i=1:nobs;
	if results.cmissing(i) == 1
	legend_colors(i,:) = hcolor(binpoints(cnt,1),:);
	cnt = cnt+1;
	else
	legend_colors(i,:) = [1 1 1];
	end;
end;
results.map_colors = legend_colors;

elseif nbars == 1
% here we have only one color, so we kludge the legend and colormap
nobs = length(variable);
legend_colors = zeros(nobs,3);
for i=1:nobs;
	if results.cmissing(i) == 1
	legend_colors(i,:) = [1 1 1];
	else
	legend_colors(i,:) = [1 1 1];
	end;
end;
results.map_colors = legend_colors;
% put a message in the legend figure window for the user
figure(results.legend_fig);
clf;
	Hwarning=uicontrol('Style','text','Units','Normalized','Position',[0.1,0.1,1,0.035],'Backgroundcolor',[1 1 1], ...
     'String',['Only 1 value for this variable'],'enable','inactive','FontSize',6,'HorizontalAlignment','left');
% warning('arc_histmap: you must select at least 2 polygons for the legend');
end;


