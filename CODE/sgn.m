function out1 = sgn(data)
%function out1 = sgn(data)
%
% the "sgn" function, which assigns the value -1 if the sign is negative, 
% 1 if the sign is positive, and 0 if the value is exactly zero.
%
%  andrew patton
%
% sat, 18 sep, 2004


out1 = (data>0) - (data<0);