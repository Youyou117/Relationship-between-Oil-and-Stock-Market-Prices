function out1 = randint(lower,upper,T,k);
%function out1 = randint(lower,upper,T,k);
%
% Generates random integers between lower and upper
%
% INPUTS:		lower, a scalar, the lower bound on draws
%					upper, a scalar, the upper bound
%					T, a scalar, the number of rows
%					k, a scalar, the number of columns
%
% OUTPUTS:	out1, a Txk matrix filled with random integers
%								in the interval [lower,upper]
%
%  Andrew Patton
%
%  15 July, 2002.

if nargin<4;
   k=1;
end
if nargin<3;
   T=1;
end
if nargin<2;
   lower=0;
   upper=1;
end
out1 = rand(T,k)*(upper-lower+1-eps)+lower-0.5+eps;	% making sure that each integer gets same prob mass
out1 = round(out1);
