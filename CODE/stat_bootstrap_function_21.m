function out1 = stat_bootstrap_function_21(T,B,length,rand_state);
%function out1 = stat_bootstrap_function_21(T,B,length);
%
% This function generates bootstrap samples of the matrix data
% and returns the time indices for each sample
%
% USING THE STATIONARY BOOTSTRAP OF POLITIS AND ROMANO (1994)
%
% INPUTS:   T, a scalar, the number of time series observations to generate
%			B, a scalar, the number of bootstrap samples
%			length, a scalar, this is the length to use for the stationary bootstrap
%           rand_state, a scalar, the "seed" for Matlab's random number generator (default is sum(100*clock) )
%
% OUTPUTS:	out1, a TxB matrix of time indices for each bootstrap sample
%
%  Andrew Patton
%
%  20 Jan 2008

if nargin<4 | isempty(rand_state)
    rand('state',sum(100*clock));		% making sure a new seed is picked each time
else
    rand('state',rand_state);
end

out1 = -999.99*ones(T,B);

% BOOTSTRAPPING THE SAMPLE AND GENERATING THE BOOTSTRAP DIST'N OF THE THREE MEASURES
for bb = 1:B;
    temp = ceil(rand*T);
    out1(1,bb) = temp;
    for tt = 2:T;
        temp2 = rand;
        if temp2>1/length		%  then we just take next obs
            if temp==T  % loop back to first obs if running past T
                temp = 1;
            else
                temp = temp+1;
            end
        else
            temp = ceil(rand*T);
        end
        out1(tt,bb) = temp;
    end
end