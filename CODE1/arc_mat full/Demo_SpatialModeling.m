
% filename = 'Shapefiles\eu138';
% 
% results = shape_read(filename);
% 
% nobs = results.nobs;
% ecu95 = results.data(:,20);
% ecu80 = results.data(:,5);
% igrwth = ecu95./ecu80;
% y = igrwth - ones(nobs,1);
%         
% latt = results.xc;
% long = results.yc;
% 
% W = make_nnw(latt,long,6); % 6 nearest neighbors weight matrix
% 
% n = length(y);
% xmat = [ones(n,1) log(ecu80) W*log(ecu80)];
% vnames = strvcat('y-growth','constant','logy80','W*logy80');
% options.vnames = vnames;
% options.nbc = 5;
% missing = ones(nobs,1);
% options.missing = missing;
% options.mapmenu = 1;
% options.legendmenu = 1;
% 
% variable = [y xmat];
% data = am_DataSource(variable,results,options);


filename = 'Shapefiles\uscounties_projected';

results = shape_read(filename);
nobs = results.nobs;
missing = ones(nobs,1);
fips = results.data(:,5);



load countyg.dat; % we load a county-level data file and re-order it
% see countyg.txt for data file documentation
id = countyg(:,1);

% order the two data sets using the fips county codes from the shape file
[nobsc,nvarsc] = size(countyg);
out = zeros(nobsc,nvarsc);
for i=1:length(fips);
fipsi = fips(i,1);
ind = find(id == fipsi);
	if length(ind) == 1
	out(i,:) = countyg(ind,:);
    else
	% fill in zeros
	out(i,:) = zeros(1,nvarsc);
	end;
end;

countyg = out;

y1 = countyg(:,4); % county employment growth rate
y2 = countyg(:,5); % county population growth rate


latt = countyg(:,2);
long = countyg(:,3);

[j W j] = xy2cont(-latt,-long);

n = length(y1);
countyg(:,16) = countyg(:,16)/1000;
xmat = [ones(n,1) countyg(:,8) countyg(:,10:end-1)];
xmat(:,2:end) = studentize(xmat(:,2:end));
vnames = strvcat('y=empgr80-90','constant','logy80','empdensity','popdensity','log area', ...
    'college','manufemp','unemploy','y-percapita','education spending','highway spending','police spending', ...
    'non-white','urban dummy');

options.vnames = vnames;
options.nbc = 5;
options.missing = missing;
options.mapmenu = 1;
options.legendmenu = 1;

variable = [y1 xmat];


data = am_DataSource(variable,results,options);

 gsar = am_GlobalSAR(data);
 lsar = am_LocalSAR(data);
 %gsar = am_GlobalSAR(data);
 %gsem = am_GlobalSEM(data);

map = am_BaseMap(data);

link =am_Linkage;
link.addMap(map);
link.addGraph(lsar);
link.addGraph(gsar);
 %link.addGraph(gsar);
% link.addGraph2(sar);
 %link.addGraph2(gsar);
 %link.addGraph(gsem);
%link.addGraph(gsar);








