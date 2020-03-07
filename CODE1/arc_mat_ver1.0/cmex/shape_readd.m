% PURPOSE: An example of using shape_read() 
%          with an ArcView shape file containing 
%          49 polygons for neighborhoods in Columbus, OH
%---------------------------------------------------
% USAGE: shape_read_d
%---------------------------------------------------

clear all;

filename = '..\shape_files\columbus';

results = shape_read(filename);

% here we show the input information
% from the shape file contained in
% the results structure variable

results

% now extract some data from results.data
% to implement a spatial regression
y = results.data(:,9);
n = results.nobs;
x = [ones(n,1) results.data(:,8) results.data(:,7)];

latt = results.yc;
long = results.xc;

[junk,W,junk] = xy2cont(latt,long);

vnames = strvcat('crime index','constant','income','house value');

sresult = sar(y,x,W);
prt(sresult,vnames);

