% PURPOSE: An example of using dbf_read() cmex function
%          reads an ArcView dbf file containing 
%          sample data for 30 provinces in China
%---------------------------------------------------
% USAGE: dbf_readd
%---------------------------------------------------

filename = '..\shape_files\michigan';

[data,vnames] = dbf_read(filename);

vnames2 = strvcat(vnames);

fprintf(1,'printing first 6 variables \n');
in.cnames = vnames2(1:6,:); % pull out first 6 variable names
in.fmt = strvcat('%12.4f');
mprint(data(1:10,1:6),in);     % pull out first 6 variable columns

fprintf(1,'size of data matrix \n');
size(data)

