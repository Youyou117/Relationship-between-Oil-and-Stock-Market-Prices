% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file containing 
%          3,142 polygons for US counties
%---------------------------------------------------
% USAGE: arc_histmapd5
%---------------------------------------------------

clear all;

filename = '../shape_files/usstates49';

results = shape_read(filename);

% results.data contains:
% col 1 AREAKEY	
% col 2 LATITUDE	
% col 3 LONGITUDE	
% col 4 POP1990	
% col 5 percapita income in 1970
% col 6 percapita income in 1980
% col 7 percapita income in 1990
% col 8 percapita income in 2000

% compute growth rates
data = results.data;
gr70to80 = log(data(:,6)) - log(data(:,5));
gr80to90 = log(data(:,7)) - log(data(:,5));
gr90to00 = log(data(:,8)) - log(data(:,7));
gr70to00 = log(data(:,8)) - log(data(:,6));

gr70to80 = gr70to80/11; % make these annualized rates
gr80to90 = gr80to90/11; % make these annualized rates
gr90to00 = gr90to00/11; % make these annualized rates
gr70to00 = gr70to00/31; % make these annualized rates

mapdata = [data(:,4:8) gr70to80 gr80to90 gr90to00 gr70to00];

vnames = strvcat('population90','pcincome70','pcincome80','pcincome90','pcincome00', ...
    'gr70to80','gr80to90','gr90to00','gr70to00');

options.vnames = vnames;

arc_histmap(mapdata,results,options);

