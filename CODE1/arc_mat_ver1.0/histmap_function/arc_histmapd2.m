% PURPOSE: An example of using arc_histmap() 
%          with an ArcView shape file containing 
%          1008 polygons for Ohio zip code areas
% demonstrates merging a data file with a shape file
%---------------------------------------------------
% USAGE: arc_histmapd2
%---------------------------------------------------

clear all;

filename = '../shape_files/ohio_zip';

results = shape_read(filename);

nobs = results.nobs;
nvar2 = results.nvars;
mzip = results.data(:,1); % we will key on this to merge 2 data sets

% variables in results.data are:
% col1 zipcode
% col2 zipcode 
% col3 latitude centroid
% col4 longitude centroid
% col5 population in 1990


% load a different data set with data on 1965 school buildings
% in 799 zip code areas
ohioschool = load('../shape_files/ohioschool.dat');
szip = ohioschool(:,1); % zip code areas in which the school buildings are located

% merge the two sets of data by zip code area
[junk,nvar1] = size(ohioschool);
nvar1 = nvar1 - 3; % we will skip the first 3 columns here because they are redundant
                   % with the information in the ohio_zip file
mdata = zeros(nobs,nvar1 + nvar2);
missing = ones(nobs,1); % an index for missing zip code areas, initialized to all non-missing
for i=1:nobs; % match ohioschool data to zip code areas in the map
zipi = mzip(i,1);
ind = find(szip == zipi);
	if length(ind) > 1
	tmp = mean(ohioschool(ind,4:end)); % where # schools > 1, use the mean
	else
	tmp = ohioschool(ind,4:end);
	end;
if length(ind) > 0
mdata(i,:) = [results.data(i,1:end) tmp];
else % if there are no schools in the zip code area use zero
mdata(i,:) = [results.data(i,1:end) zeros(1,nvar1)];
missing(i,1) = 0; % make this zip code area a missing value
end;
end;
        


% variables in ohioschools.data are:
% data for 2001-02 year
% col 1 = zip code
% col 2 = lattitude (zip centroid)
% col 3 = longitude (zip centroid)
% col 4 = buidling irn
% col 5 = district irn
% col 6 = # of teachers (FTE 2001-02)
% col 7 = teacher attendance rate
% col 8 = avg years of teaching experience
% col 9 = avg teacher salary
% col 10 = Per Pupil Spending on Instruction
% col 11 = Per Pupil Spending on Building Operations
% col 12 = Per Pupil Spending on Administration
% col 13 = Per Pupil Spending on Pupil Support
% col 14 = Per Pupil Spending on Staff Support
% col 15 = Total Expenditures Per Pupil
% col 16 = Per Pupil Spending on Instruction % of Total Spending Per Pupil
% col 17 = Per Pupil Spending on Building Operations % of Total Spending Per Pupil
% col 18 = Per Pupil Spending on Administration % of Total Spending Per Pupil
% col 19 = Per Pupil Spending on Pupil Support % of Total Spending Per Pupil
% col 20 = Per Pupil Spending on Staff Support % of Total Spending Per Pupil
% col 21 = irn number
% col 22 = avg of all 4th grade proficiency scores
% col 23 = median of 4th grade prof scores
% col 24 = building enrollment
% col 25 = short-term students < 6 months
% col 26 = 4th Grade (or 9th grade) Citizenship % Passed 2001-2002
% col 27 = 4th Grade (or 9th grade)  math % Passed 2001-2002
% col 28 = 4th Grade (or 9th grade)  reading % Passed 2001-2002
% col 29 = 4th Grade (or 9th grade)  writing % Passed 2001-2002
% col 30 = 4th Grade (or 9th grade)  science % Passed 2001-2002
% col 31 = pincome per capita income in the zip code area
% col 32 = nonwhite percent of population that is non-white
% col 33 = poverty percent of population in poverty
% col 34 = samehouse % percent of population living in same house 5 years ago
% col 35 = public % of population attending public schools
% col 36 = highschool
% col 37 = assoc educ attainment for 25 years plus
% col 38 = college
% col 39 = grad
% col 40 = prof

scores = mdata(:,29); % col 22 plus 7 columns from results.data = col 29


vnames = strvcat('zipcode','zipcode','latt','long','pop90', ...
'buildingirn','districtirn','#teachers','teacher attendance','teacher experience','teacher salary', ...
'per pupil spending instruction','per pupil spending buidling','per pupil speinding on administration',...
'per pupil spending pupil support', ...
'per pupil spending staff support','total expend per pupil','% on instruction','% on building',...
'% on administration','% on pupil support','% on staff support','irn number','avg 4th or 9th grade scores',...
'median 4th or 9th grade scores','building enrollment','short-term students','4th or 9th grade citizenship',...
'4th or 9th grade math','4th or 9th grade reading','4th or 9th grade writing','4th or 9th grade science','per captia income',...
'nonwhite','poverty','samehouse','public','highschool','associate degree','college','graduate','professional');


options.vnames = vnames;
options.nbc = 5;
options.missing = missing;
options.mapmenu = 1;
options.legendmenu = 1;
arc_histmap(mdata,results,options);

