% merge projected county shape file info 
clear all;

filename = '..\shape_files\us_counties';

tic;
results = shape_read(filename);
toc;

% results.data contains:
% col 1 blank area name
% col 2 fips code	
% col 3 LATITUDE	
% col 4 LONGITUDE	
% col 5 POP1990	
% col 6 population living in same house 5 years ago
% col 7 POP65PLUS
% col 8 EMPLOYMENT	
% col 9 HOUSEHLDS	
% col 10 MEDINCOME	
% col 11 PCINCOME	
% col 12 POVERTY	
% col 13 NOHOUSES	
% col 14 MEDRENT	
% col 15 WHITEPOP	
% col 16 BLACKPOP	
% col 17 NATIVPOP

% pull out some variables to map
mapdata = results.data(:,3:end);
% mapdata matrix contains
% col 1 LATITUDE	
% col 2 LONGITUDE	
% col 3 POP1990	
% col 4 population living in same house 5 years ago
% col 5 POP65PLUS
% col 6 EMPLOYMENT	
% col 7 HOUSEHLDS	
% col 8 MEDINCOME	
% col 9 PCINCOME	
% col 10 POVERTY	
% col 11 NOHOUSES	
% col 12 MEDRENT	
% col 13 WHITEPOP	
% col 14 BLACKPOP	
% col 15 NATIVPOP


% convert whitepop to a percentage of all pop
mapdata(:,13) = mapdata(:,13)./mapdata(:,3);
% convert blackpop to a percent of all pop
mapdata(:,14) = mapdata(:,14)./mapdata(:,3);
% convert native american pop to a percent of all pop
mapdata(:,15) = mapdata(:,15)./mapdata(:,3);
% convert population 65 plus to a percentage of all pop
mapdata(:,5) = mapdata(:,5)./mapdata(:,3);
% convert employment to be a percentage of all pop
mapdata(:,6) = mapdata(:,6)./mapdata(:,3);
% convert those living in the same house 5 years ago to be a percent of all houses
mapdata(:,4) = mapdata(:,4)./mapdata(:,11);

% pull out associated variable names
vnames = results.vnames(3:end,:);

options.vnames = vnames;
options.legendmenu = 1; % place a menu on the map legend for printing/editing/saving purposes
options.mapmenu = 1; % place a menu on the map for printing/editing/saving purposes

arc_histmap(mapdata,results,options);

