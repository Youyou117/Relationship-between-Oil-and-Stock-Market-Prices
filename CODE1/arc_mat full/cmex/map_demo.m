% PURPOSE: An example of using shape_read() and make_map()
%          reads an ArcView shape file, then creates a map
% NOTE: the make_map function is used by: arc_histmap, arc_moranmap, arc_sarmap
% see arc_histmapd, arc_moranmapd, arc_sarmapd for examples
%     of using these NICE mapping functions 
% This file is intended as a demo for hardcore programmers, not a typical user.
%---------------------------------------------------
% USAGE: map_demo
%---------------------------------------------------

clear all;

filename = '..\shape_files\columbus';

results = shape_read(filename);
poly = make_map(results);

% add colors to each polygon based on the variable crime
crime = results.data(:,9);

hc = colormap('hsv');
% hc is always 64 by 3 matrix

nbc = 5;
% hc is always 64 by 3 matrix
incr = floor(64/nbc);
cindex = 1:incr:64;
cindex = cindex(1:nbc);
hcolor = hc(cindex,:);

grid = (max(crime)-min(crime))/nbc;

c1 = find(crime <= grid);
c2 = find(crime > grid & crime <= 2*grid);
c3 = find(crime > 2*grid & crime <= 3*grid);
c4 = find(crime > 3*grid & crime <= 4*grid);
c5 = find(crime > 4*grid);
chk = [c1' c2' c3' c4' c5'];

for i=1:results.npoly;
if length(find(c1 == i)) > 0
	for k=1:results.nparts(i);
	set(poly(i).handles(k),'FaceColor',hcolor(1,:));
	end;
elseif length(find(c2 == i)) > 0
	for k=1:results.nparts(i);
    set(poly(i).handles(k),'FaceColor',hcolor(2,:));
	end;
elseif length(find(c3 == i)) > 0
	for k=1:results.nparts(i);
    set(poly(i).handles(k),'FaceColor',hcolor(3,:));
	end;
elseif length(find(c4 == i)) > 0
	for k=1:results.nparts(i);
    set(poly(i).handles(k),'FaceColor',hcolor(4,:));
	end;
elseif length(find(c5 == i)) > 0
	for k=1:results.nparts(i);
    set(poly(i).handles(k),'FaceColor',hcolor(5,:));
	end;
end;
end;

% turn on the map using the poly(1).fig_handle
set(poly(1).fig_handle,'Visible','on');
pause;

clf reset;
close all;