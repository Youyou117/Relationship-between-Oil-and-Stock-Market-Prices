function [Gcdf, Gpdf, x_grid, sig_grid] = Gcdfpdf_wholeY_Skewtt(q,theta,Weight_GL_Ginv,Npoints)
sig_grid = [0.001,0.01,(0.1:0.1:2),(2.1:0.1:5.6)]';

x1 = -15;
x2 = 15;

x_grid = nan(Npoints,1);

if length(x_grid(:,1)) ~= length((x1:(x2-x1)/(Npoints-1):x2)')
    x_grid(:,1) = [(x1:(x2-x1)/(Npoints-1):x2)';x2];
else
    x_grid(:,1) = (x1:(x2-x1)/(Npoints-1):x2)';
end

x_grid = [-30;x_grid;30];

Gcdf = nan(length(x_grid),length(sig_grid));
Gpdf = nan(length(x_grid),length(sig_grid));

for ii=1:length(x_grid);
    for jj = 1:length(sig_grid);
    Gpdf(ii,jj) = factor_cop_Gpdf_GLQuad_Skewtt(x_grid(ii),[sig_grid(jj);theta(2);theta(3);theta(4)],Weight_GL_Ginv);
    Gcdf(ii,jj) = factor_cop_Gcdf_GL_Skewtt(x_grid(ii),[sig_grid(jj);theta(2);theta(3);theta(4)],Weight_GL_Ginv);
    end
end
