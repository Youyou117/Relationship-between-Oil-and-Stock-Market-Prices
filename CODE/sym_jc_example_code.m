% Some sample code showing how I use my Symmetrised Joe-Clayton
% copula programs
%
%  Andrew Patton
%
%  Wednesday, 12 February, 2003.

tauU = 0.2;
tauL = 0.6;
T = 1000;

tic;
data = sym_jc_rnd(tauU,tauL,T);		% generating some data from the symmetrised Joe-Clayton copula
toc	% takes 1.7 seconds for T=1000 on my computer (2.4Ghz)
		% takes 7.5 seconds for T=5000 on my computer (2.4Ghz)

%figure(1),scatter(data(:,1),data(:,2))		% plotting the data
%figure(2),scatter(norminv(data(:,1)),norminv(data(:,2)))		% plotting the data transformed to have N(0,1) margins
% evidence of asymmetry is easier to see in second graph than first (at least to my eyes)


% Maximum likelihood estimation
lower = 1e-12*ones(2,1);			% setting lower bound just above 0
upper = 1-1e-12*ones(2,1);			% upper bound just below 1
theta0 = [0.4;0.4];					% starting values
thetahat  = fmincon('sym_jc_cl',theta0,[],[],[],[],lower,upper,[],[],data); 
[[tauU;tauL],thetahat]
thetahat2 = fmincon('sym_jc_cl',[0.24;0.75],[],[],[],[],lower,upper,[],[],data); 
[[tauU;tauL],thetahat2]



