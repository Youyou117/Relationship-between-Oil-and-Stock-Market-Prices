% Example code to illustrate main function for running "bounds" tests
%
% Used in the papers:
%
% Patton, A.J., and A. Timmermann, 2012, Forecast Rationality Tests Based on Multi-Horizon Bounds, 
% Journal of Business and Economic Statistics, 30(1), 1-17.
% Paper available at: http://public.econ.duke.edu/~ap172/Patton_Timmermann_bounds_JBES_2012.pdf
%
%  Andrew Patton
%
%  22 December 2010


% First step - compute some simulated multi-horizon forecasts. I will use an AR(1)
muy = 0.75;
phi1 = 0.5;
sig2y = 0.5;
sig2eps = sig2y*(1-(phi1^2));
T = 100;  % sample size
Hmax = 8;  % maximum forecast horizon

% simulating the target variable
y = nan(T,1);
y(1) = muy;
for tt=2:T;
    y(tt) = muy + phi1*(y(tt-1)-muy) + randn*sqrt(sig2eps);
end
figure(1),plot(y);
figure(2),sacf(y,10);

% generating the forecasts, lined up in "event time", so that row t containts yhat[t|t-h]
yhat = nan(T,Hmax);
for hh=1:Hmax;
    yhat(1+hh:end,hh) = muy + (phi1.^hh).*(y(1:end-hh)-muy);
end
figure(3),plot(yhat),legend(int2str((1:Hmax)'))

% plotting the MSE and forecat variance
MSE = nanmean(  (y*ones(1,Hmax)-yhat).^2 )';
MSF = diag(nancov(  yhat ));
figure(4),bar(-(Hmax:-1:1)',[MSE(end:-1:1),MSF(end:-1:1)],'stacked'),colormap('summer');hold on;
plot(-(Hmax+1:-1:0)',cov(y)*ones(Hmax+2,1),'k--','LineWidth',2),...
    axis([-Hmax-0.5,-0.5,0,1.33*cov(y)]),legend('MSE','V[yhat]','V[y]'),xlabel('Forecast horizon'),title('MSE and MSF for a simulated AR(1)');hold off;


% now running the test (results will be printed to screen in a little table)
tic;[out1, out2] = bounds_tests_function_40nan_simple(y,yhat,1); toc  % takes about 0.8 seconds on my 2012 Lenovo laptop.