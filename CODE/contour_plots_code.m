% Code to generate iso-probability contour plots of bivariate
% distributions, all with N(0,1) marginal distributions and connected via various copulas 
%
%  Andrew Patton
%
%  21 august 2006

% Written for the following paper:
%
% Patton, A.J., 2006, Copula-Based Models for Financial Time Series. 
% To be published in T.G. Andersen, R.A. Davis, J.-P. Kreiss and T. Mikosch (eds.), 
% "Handbook of Financial Time Series", Springer Verlag. 
%
% http://fmg.lse.ac.uk/~patton



T = 100;  % takes only around 15 seconds when T=100 on a single-processor 2.6GHz machine

% i cannot see a difference in the printed graph based on T=100 vs T=500,
% though it is a lot smaller and faster to load. so just use T=100;

tic;
xx = (-2:4/(T-1):2)';
uu = normcdf(xx);

v = (0.02:0.03:0.2);

% 1. Normal copula
rho = 0.5;
zz = nines(T,T);  % this the the part of the pdf from the X variable
tic;
for ii=1:T;
   zz(:,ii) = bivnormpdf(xx,xx(ii),[0,0],theta2rho(rho));
end
toc

figure(101);subplot(3,2,1),contour(xx,xx,zz,v,'r-'),...
    title('Normal copula, \rho = 0.5');
figure(102);subplot(3,2,1),contour(xx,xx,zz,v,'k-'),...
    title('Normal copula, \rho = 0.5');
hold on;subplot(3,2,1),plot([-2,2],[0,0],'k--',[0,0],[-2,2],'k--');hold off;    

% 2. Student's t copula
rho = 0.5;
nu = 3;
zz = normpdf(xx)*ones(1,T);  % this the the part of the pdf from the X variable
tic;
for ii=1:T;
   zz(:,ii) = zz(:,ii).*normpdf(xx(ii)).*tcopula_pdf(uu,uu(ii),rho,nu);
end
toc
figure(101);subplot(3,2,2),contour(xx,xx,zz,v,'m-'),...
    title('Student''s t copula, \rho = 0.5, \nu = 3'),...
figure(102);subplot(3,2,2),contour(xx,xx,zz,v,'k-'),...
    title('Student''s t copula, \rho = 0.5, \nu = 3'),...
hold on;subplot(3,2,2),plot([-2,2],[0,0],'k--',[0,0],[-2,2],'k--');hold off;    

% 3. Clayton copula
kappa = 1;
zz = normpdf(xx)*ones(1,T);  % this the the part of the pdf from the X variable
tic;
for ii=1:T;
   zz(:,ii) = zz(:,ii).*normpdf(xx(ii)).*clayton_pdf(uu,uu(ii),kappa);
end
toc

figure(101);subplot(3,2,3),contour(xx,xx,zz,v,'Color',[0 0.7 0]),...
    title('Clayton copula, \kappa = 1');
figure(102);subplot(3,2,3),contour(xx,xx,zz,v,'k-'),...
    title('Clayton copula, \kappa = 1');
hold on;subplot(3,2,3),plot([-2,2],[0,0],'k--',[0,0],[-2,2],'k--');hold off;    


% 4. Gumbel copula
kappa = 1.5;
zz = normpdf(xx)*ones(1,T);  % this the the part of the pdf from the X variable
tic;
for ii=1:T;
   zz(:,ii) = zz(:,ii).*normpdf(xx(ii)).*gumbel_pdf(uu,uu(ii),kappa);
end
toc

figure(101);subplot(3,2,4),contour(xx,xx,zz,v,'Color',[1 0.4 0]),...
    title('Gumbel copula, \kappa = 1.5');
figure(102);subplot(3,2,4),contour(xx,xx,zz,v,'k-'),...
    title('Gumbel copula, \kappa = 1.5');
hold on;subplot(3,2,4),plot([-2,2],[0,0],'k--',[0,0],[-2,2],'k--');hold off;    


% 5. SJC copula
tauU = 0.45;
tauL = 0.2;
zz = normpdf(xx)*ones(1,T);  % this the the part of the pdf from the X variable
tic;
for ii=1:T;
   zz(:,ii) = zz(:,ii).*normpdf(xx(ii)).*sym_jc_pdf(uu,uu(ii),tauU,tauL);
end
toc
figure(101);subplot(3,2,5),contour(xx,xx,zz,v,'Color',[0.5 0.5 0.5 ]),...
    title('SJC copula, \tau^U = 0.45, \tau^L = 0.2');
figure(102);subplot(3,2,5),contour(xx,xx,zz,v,'k-'),...
    title('SJC copula, \tau^U = 0.45, \tau^L = 0.2');
hold on;subplot(3,2,5),plot([-2,2],[0,0],'k--',[0,0],[-2,2],'k--');hold off;    


% 6. Mixture of normals copula
rho1 = 0.95;
rho2 = 0.05;
zz = nines(T,T);  % this the the part of the pdf from the X variable
tic;
for ii=1:T;
   zz(:,ii) = 0.5*bivnormpdf(xx,xx(ii),[0,0],theta2rho(rho1)) + 0.5*bivnormpdf(xx,xx(ii),[0,0],theta2rho(rho2));
end
toc
figure(101);subplot(3,2,6),contour(xx,xx,zz,v,'Color',[0 0 1 ]),...
    title('Mixed normal copula, \rho_1 = 0.95, \rho_2 = 0.05');
figure(102);subplot(3,2,6),contour(xx,xx,zz,v,'k-'),...
    title('Mixed normal copula, \rho_1 = 0.95, \rho_2 = 0.05');
hold on;subplot(3,2,6),plot([-2,2],[0,0],'k--',[0,0],[-2,2],'k--');hold off;    

