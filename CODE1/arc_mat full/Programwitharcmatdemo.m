filename = 'Shapefiles\eu138';

results = shape_read(filename);

nobs = results.nobs;
ecu95 = results.data(:,20);
ecu80 = results.data(:,5);
igrwth = ecu95./ecu80;
y = igrwth - ones(nobs,1);
        
latt = results.xc;
long = results.yc;

W = make_nnw(latt,long,6); % 6 nearest neighbors weight matrix

n = length(y);
xmat = [ones(n,1) log(ecu80) W*log(ecu80)];
vnames = strvcat('y-growth','constant','logy80','W*logy80');
options.vnames = vnames;
options.nbc = 5;
missing = ones(nobs,1);
options.missing = missing;
options.mapmenu = 1;
options.legendmenu = 1;

variable = [y xmat];


data = am_DataSource(variable,results,options);


%data = am_DataSource(variable,results,options); %create data source

map = am_BaseMap(data);%create maps and statistical plots
moran = am_Moranplot(data);
gsar = am_GlobalSAR(data);
%m3d = am_Scatter3d(data);
hist = am_Histogram(data);

link = am_Linkage;%link all plots together
link.addMap(map);
link.addGraph(hist);
link.addGraph(moran);
%link.addGraph2(m3d);
link.addGraph(gsar);
