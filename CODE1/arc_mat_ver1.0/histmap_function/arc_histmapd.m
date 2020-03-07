% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file containing 
%          30 polygons for China provinces 
%---------------------------------------------------
% USAGE: arc_histmapd
%---------------------------------------------------

clear all;

filename = '../shape_files/china';
% see china.txt for data file documentation

results = shape_read(filename);

% we pull out selected variables and form a matrix of things to map
gdp80 = 1000*results.data(:,17); % 1980 GDP
gdp95 = 1000*results.data(:,20); % 1995 GDP
pop95 = results.data(:,5); % 1995 population
pop80 = results.data(:,2); % 1980 population

gdppop80 = gdp80./pop80; % per capita gdp 1980
gdppop95 = gdp95./pop95; % per capita gdp 1995

gdpgrwth = log(gdppop95) - log(gdppop80); % growth of gdp per capita 1980-1995
gdpgrowth = gdpgrwth/15;
gdpdiff = gdppop95 - gdppop80;            % change in gdp per capita 1980-1995

% our matrix of variables to map
pltvariables = [gdpgrwth gdppop80 gdppop95 gdpdiff];

options.vnames = strvcat('gdp/pop growth 80-95','gdp/pop80','gdp/pop95','gdp/pop difference 80-95');

arc_histmap(pltvariables,results,options);

