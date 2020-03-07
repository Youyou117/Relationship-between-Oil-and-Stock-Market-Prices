% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file that was saved 
%          with Albers Equal-Area Conic projection
%          (a typical US map projection)
% see the ArcView (version 3.2a) project file:
% ..\shapefiles\uscounties_projected.apr
%---------------------------------------------------
% USAGE: arc_histmapd9
%---------------------------------------------------

clear all;

filename = '../shape_files/uscounties_projected';

results = shape_read(filename);

% results.data contains:
% col1 AREANAME
% col2 FIPS
% col3 LATITUDE
% col4 LONGITUDE
% col5 POP1990
% col6 1987_PCI
% col7 1988_PCI
% col8 1989_PCI
% col9 1990_PCI
% col10 1991_PCI
% col11 1992_PCI
% col12 1993_PCI

% clever way to produce a state index
fips = results.data(:,2);
states = round(fips/1000);

% we add the state index 
% and pass everything along to the mapping function
data = [states results.data];
options.vnames = strvcat('states',results.vnames);

arc_histmap(data,results,options);

