function results = am_moran(results)
% PURPOSE: compute moran scatterplot figure
% returns an n x 3 matrix of colors that can be used to label the map polygons
% and a handle to the legend figure

cmap = results.cmap;
nbc = results.nbc;
vindex = results.vindex;
vflag = results.vflag;

cvariable = results.cvariable; 
WX = results.WX(:,vindex);
obs_selected = results.obs_selected;
nobs = length(obs_selected);
if vflag == 1
vname = results.vnames(vindex,:);
elseif vflag == 0
vname = results.vnames(vindex+1,:);
end;

Q0 = results.Q0;
Q1 = results.Q1;
Q2 = results.Q2;
Q3 = results.Q3;
Q4 = results.Q4;


if (results.moran_fig == 0)
	if results.legendmenu == 0
	moran_fig = figure('Position',[results.width+60 100 400 400], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Moran Scatterplot', ...
                      'MenuBar','none');
	elseif results.legendmenu == 1
	moran_fig = figure('Position',[results.width+60 100 400 400], ... % [left bottom width height]
                      'NumberTitle','off', ...
                      'Name','Moran Scatterplot');
	end;
    


results.moran_fig = moran_fig;
else
figure(results.moran_fig);
end;

hc = colormap(cmap);
% hc is always 64 by 3 matrix
incr = floor(64/nbc);
cindex = 1:incr:64;
cindex = cindex(1:nbc);
hcolor = hc(cindex,:);

map_colors = ones(nobs,3);
for i=1:3;
map_colors(Q0,i) = 1;
map_colors(Q1,i) = hcolor(1,i);
map_colors(Q2,i) = hcolor(2,i);
map_colors(Q3,i) = hcolor(3,i);
map_colors(Q4,i) = hcolor(4,i);
end;
results.map_colors = map_colors;

if nobs > 1000
msize = 4;
elseif nobs > 100
msize = 8;
elseif nobs > 50
msize = 16;
else
msize = 32;
end;

good = find(results.cmissing == 1);

mycolor = zeros(nobs,3);
mycolor(Q0,:) = ones(length(Q0),3);
mycolor(Q1,:) = matmul(hcolor(1,:),ones(length(Q1),3));
mycolor(Q2,:) = matmul(hcolor(2,:),ones(length(Q2),3));
mycolor(Q3,:) = matmul(hcolor(3,:),ones(length(Q3),3));
mycolor(Q4,:) = matmul(hcolor(4,:),ones(length(Q4),3));

hscattergroup = scatter(cvariable,WX,msize,mycolor); %get the handle of the scattergroup object created

hpatch = get(hscattergroup,'Children');  %return the handles (array m_fig) of patches. These patchs i.e dots on the scatterplot
                                        %are scattergroup' children.

% after test, I found that the ith handle corresponds to the jth region on the map, and the
% relationship between i and j is i+j = the total number of regions + 1

set(hpatch(Q0),'Visible','off'); %hide all dots
set(hpatch(Q1),'Visible','off');
set(hpatch(Q2),'Visible','off');
set(hpatch(Q3),'Visible','off');
set(hpatch(Q4),'Visible','off');

axis normal;

%trans help to find out the right handles according to obs_selected
trans = ones(length(obs_selected),1);
trans = trans*(length(Q0)+length(Q1)+length(Q2)+length(Q3)+length(Q4)+1);
trans = matsub(trans,obs_selected);

set(hpatch(trans),'Visible','on');

set(gcf,'UserData',hpatch);

% try keying this only on selected observations
% for the case of a zoom
tmpvariable = cvariable(results.obs_selected,1);
tmpWX = WX(results.obs_selected,1);

xbuffer = 0.5*abs(0.1*min(tmpvariable));
ybuffer = 0.1*min(tmpWX);

maxx = max(tmpvariable);
minx = min(tmpvariable);
maxy = max(tmpWX);
miny = min(tmpWX);

axis([1.1* minx 1.1*maxx 1.1*miny 1.1*maxy]);

line([1.1*minx 1.1*maxx],[0 0],'Color','black');
line([0 0],[1.1*miny 1.1*maxy],'Color','black');
axis manual;

  if results.labels == 1
	for i=1:length(obs_selected);
    j = obs_selected(i,1);
	hi = text(cvariable(j)+xbuffer,WX(j),num2str(j));
	set(hi,'fontsize',6);
	end;
  end;

set(results.moran_fig,'Visible','on');
xlabel(vname);
mname = ['W*',vname];
ylabel(mname);


