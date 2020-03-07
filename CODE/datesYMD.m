function out1 = datesYMD(date)
% function out1 = datesYMD(date);
%
% INPUT:	date, a Tx1 vector of dates in the form YYYYMMDD
%
% OUTPUT:	out1, a Tx3 matrix of dates = [YYYY,MM,DD]
%
% Andrew Patton, 2002.

out1      = nan(size(date,1),3);
out1(:,1) = (date-mod(date,10000))/10000;
out1(:,2) = (mod(date,10000)-mod(date,100))/100;
out1(:,3) = mod(date,100);
