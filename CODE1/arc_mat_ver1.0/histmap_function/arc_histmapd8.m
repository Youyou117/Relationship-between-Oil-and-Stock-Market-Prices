% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file containing 
%          68 polygons for counties in Michigan
%---------------------------------------------------
% USAGE: arc_histmapd8
%---------------------------------------------------

clear all;

filename = '../shape_files/michigan';

results = shape_read(filename);

% we just pass everything along to the mapping function
data = results.data;
options.vnames = results.vnames;

arc_histmap(data,results,options);

