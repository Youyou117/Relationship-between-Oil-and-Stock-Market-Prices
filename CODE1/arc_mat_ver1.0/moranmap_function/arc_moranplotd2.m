% PURPOSE: An example of using arc_moranplot() 
%          with an ArcView shape file containing 
%          3,111 polygons for US counties
%---------------------------------------------------
% USAGE: arc_moranplotd
%---------------------------------------------------

clear all;

filename = 'uscounties_projected';

results = shape_read(filename);

% col 1  AREANAME 
% col 2  FIPS     
% col 3  LATITUDE 
% col 4  LONGITUDE
% col 5  POP1990  
% col 6  1987_PCI (per capita income)
% col 7  1988_PCI 
% col 8  1989_PCI 
% col 9  1990_PCI 
% col 10 1991_PCI 
% col 11 1992_PCI 
% col 12 1993_PCI 


    
vnames = strvcat('pop1990','1987 PCI','1988 PCI','1989 PCI','1990 PCI','1991 PCI', ...
'1992 PCI','1993 PCI');

map_data = results.data(:,5:end);

map_data(:,1) = log(map_data(:,1));

latt = results.data(:,3);
long = results.data(:,4);

% make_nnw is a function from the econometric toolbox
% that constructs a spatial weight matrix based on the
% 5 nearest neighbors
W = make_neighborsw(latt,long,6); 


% we pass everything along to the mapping function
options.vnames = vnames;

arc_moranplot(map_data,W,results,options);

