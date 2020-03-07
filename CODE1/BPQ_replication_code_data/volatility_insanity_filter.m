function [out1,count] = volatility_insanity_filter(forecasts,LB,UB,adj_forecast)
% function [out1,count] = volatility_insanity_filter(forecasts,LB,UB,adj_forecast)
%
% Function to take in a set of volatility forecasts, and change any
% forecasts that are below LB or above UB to simply adj_forecast
%
%  INPUTS:  forecasts, a Tx1 vector of forecasts
%           LB, a scalar or a Tx1 vector, the lower bound for a "reasonable" forecast. (Set to -Inf to impose no lower bound.)
%           UB, a scalar or a Tx1 vector, the upper bound for a "reasonable" forecast. (Set to +Inf to impose no lower bound.)
%           adj_forecast, a scalar, the value to use for all "unreasonable" forecasts
%
%  OUTPUTS: out1, a Tx1 vector of adjusted forecasts
%           count, a scalar, the number of forecasts that were adjusted
%
%  Andrew Patton
%
%  18 Oct 2016


out1 = forecasts;
count = sum( (forecasts<LB) + (forecasts>UB) );
out1(forecasts<LB)=adj_forecast;
out1(forecasts>UB)=adj_forecast;
