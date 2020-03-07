% compute growth rate of
% University of Toledo 
% enrollment by states and 88 Ohio counties

clear all;

mresults = shape_read('..\shape_files\usstates49');
mfips = mresults.data(:,1)/1000;

[a,b] = xlsread('enroll_by_state.xls');
n = 49; % 48 states plus DC

sfips = a(2:n+1,1);
enroll98 = a(2:n+1,2);
enroll03 = a(2:n+1,3);

% arrange this data the same as the map data
out = zeros(n,2);
for i=1:n;
    fipi = sfips(i,1);
    ind = find(mfips == fipi);
    if length(ind) > 0
        out(ind,:) = [enroll98(i,1) enroll03(i,1)];
    end;
end;

enroll98 = out(:,1);
enroll03 = out(:,2);

egrowth = log(enroll03) - log(enroll98);
egrowth = egrowth/5;


W = make_neighborsw(mresults.xc,mresults.yc,5);

% put together a matrix of information for mapping
map = [egrowth*100 log(enroll98) log(enroll03)];
options.vnames = strvcat('enroll growth 98-05','log(enroll98)','log(enroll03)');

arc_histmap(map,mresults,options);
pause;

% do moran scatterplot map of map matrix
options.labels = 1;
arc_moranplot(map,W,mresults,options);




