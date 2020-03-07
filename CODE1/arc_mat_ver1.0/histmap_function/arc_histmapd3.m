% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file containing 
%          138 polygons for European union regions
%---------------------------------------------------
% USAGE: arc_histmapd8
%---------------------------------------------------

clear all;

filename = '../shape_files/eu138';

results = shape_read(filename);

nobs = results.nobs;
% pull out some variables to map
ecu95 = results.data(:,20);
ecu80 = results.data(:,5);
igrwth = ecu95./ecu80;
igrwth = igrwth - ones(nobs,1);        
igrwth = igrwth/16; % annual percentage growth in per capita income

variables = [igrwth ecu80 ecu95]; % our mapping variables matrix
vnames = strvcat('pcapita income growth','ecu80','ecu95');
options.vnames = vnames;

arc_histmap(variables,results,options);

