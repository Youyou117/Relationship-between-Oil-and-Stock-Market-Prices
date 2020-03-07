% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file containing 
%          2862 polygons for census tracts in Ohio
%---------------------------------------------------
% USAGE: arc_histmapd7
%---------------------------------------------------

clear all;

filename = '../shape_files/ohio_tracts';

tic;
results = shape_read(filename);
toc;

% variables in results.data are:
% col 1 AREANAME	
% col 2 AREAKEY	
% col 3 LATITUDE	
% col 5 LONGITUDE	
% col 6 POP1990	
% col 7 MEDINCOME	
% col 8 PERCAPINCO	
% col 9 POVERTY
        
% make all values below 40 degrees latitude missing values
ind = find(results.data(:,3) > 40.0);
missing = zeros(results.nobs,1);
missing(ind,1) = 1;

options.vnames = results.vnames;
options.missing = missing;
mapdata = results.data;

arc_histmap(mapdata,results,options);
