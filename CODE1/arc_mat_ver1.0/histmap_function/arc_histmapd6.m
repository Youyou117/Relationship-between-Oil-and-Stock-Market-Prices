% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file containing 
%          49 polygons for neighborhoods in Columbus, OH
%---------------------------------------------------
% USAGE: arc_histmapd6
%---------------------------------------------------

clear all;

filename = '../shape_files/columbus';

results = shape_read(filename);

% here we demonstrate just feeding down the input information
% from the shape file
options.vnames = results.vnames;

arc_histmap(results.data,results,options);

