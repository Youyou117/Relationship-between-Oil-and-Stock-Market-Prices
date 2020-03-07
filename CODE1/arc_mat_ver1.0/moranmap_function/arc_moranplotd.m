% PURPOSE: An example of using arc_moranplot() 
%          with an ArcView shape file containing 
%          30 polygons for China provinces 
%---------------------------------------------------
% USAGE: arc_moranplotd
%---------------------------------------------------

%clear all;

filename = '../shape_files/china';
% see china.txt for data file documentation

results = shape_read(filename);

nobs = results.nobs;
pop95 = results.data(:,5); % 1995 population
pop80 = results.data(:,2); % 1980 population
pgrwth = pop95./pop80;
pgrwth = pgrwth - ones(nobs,1); % percentage growth over the 80 to 95 period

variable = [pop95 pop80 pgrwth];

latt = results.xc;
long = results.yc;

[j,W,j] = xy2cont(long,latt);

options.vnames = strvcat('pop95','pop80','pop growth');
%options.labels = 1;
%options.mapmenu = 1;
%options.legendmenu = 1;

arc_moranplot(variable,W,results,options);

