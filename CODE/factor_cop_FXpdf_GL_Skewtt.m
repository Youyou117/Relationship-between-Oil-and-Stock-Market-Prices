function out1 = factor_cop_FXpdf_GL_Skewtt(x, theta, GLweight, group_code, epsi)
% function out1 = factor_cop_FXpdf_GL_Skewtt(x, theta, GLweight, group_code, epsi)
%
% The joint density of [X1,...,XN] associated with skew t - t factor model (1st element of out1)
% with other joint densities evaluated at (factor loading + [0,..,eps_i,..,0], nuinv_z, nuinv_esp, psi_z)
%
% INPUTS: x,    a [N by 2] matrix, the value of x that we will evaluate g(x1,...,xN) at
%               1st column: Ginv(u) evaluated at theta (= [factor loading, nuinv_z, nuinv_esp, psi_z])
%               2nd column: Ginv(u) evaluated at theta (= [(factor loading + [0,..,eps_i,..,0]), nuinv_z, nuinv_esp, psi_z]) 
%         theta,      =[factor loading, nuinv_z, nuinv_esp, psi_z], the parameters of Fz and Feps
%         GLweight,   a matrix, nodes and weight for Gauss-Legendre quadrature 
%         group_code, a Nx1 vector of group codes into which each firm is classified 
%         epsi,       a scalar, the step size for numerical derivative for score calculation
%
%  OUTPUTS: out1, a  [(1 + # of group) by 1] vector;
%                 1st element: the value of the joint pdf of skew t-t evaluated at
%                  theta (= [factor loading, nuinv_z, nuinv_esp, psi_z])
%                 (i+1)th elements: the value of the joint pdf of skew t-t evaluated at
%                  theta (= [factor loading + [0,..,eps_i,..,0], nuinv_z, nuinv_esp, psi_z])
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 April 2015
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence 
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University. 


out1 = GLquad(@factor_cop_FXpdf_calc1_Skewtt_DiffLoad_VT, 1e-5, 1-1e-5, GLweight, x, theta, group_code, epsi);


    function out = factor_cop_FXpdf_calc1_Skewtt_DiffLoad_VT(u, x, theta, group_code, epsi)
        
        % out = Nnodes * (1+Ngroup)
        
        Ngroup = max(group_code);
        if length(theta(1:end-3)) ~= Ngroup
            error('N_group is not equal to N_theta');
        end
        
        nuinv_z     = theta(end-2);
        nuinv_eps   = theta(end-1);
        psi_z       = theta(end);
        
        N = size(x,1);
        [Nnodes,k] = size(u);  % skewtdis_inv needs "u" to be a column vector rather than a row vector
        if k>Nnodes
            u=u';
            [Nnodes,~] = size(u);
        end
                
        Fz_inv_u = skewtdis_inv(u, 1/nuinv_z, psi_z);
        
        xxx = nan(Nnodes*N,1);
        for ii=1:N;
            inx = group_code(ii);
            xxx( Nnodes*(ii-1)+1 : Nnodes*ii, 1)= x(ii,1) - sqrt(theta(inx))*Fz_inv_u; % Nnodes * 1
        end
        
        xxx = xxx(:, ones(1,Ngroup+1)); 
        
        for ii=1:N;
            inx = group_code(ii) ;
            xxx( Nnodes*(ii-1)+1 : Nnodes*ii, 1+inx) = x(ii,2) - sqrt(theta(inx)+epsi)*Fz_inv_u;
        end
        
        xxx = xxx(:);
        
        out_temp = skewtdis_pdf( xxx, 1/nuinv_eps, 0 );
        out_temp = reshape(out_temp, Nnodes*N,(Ngroup+1));
                
        out = nan(Nnodes,Ngroup + 1);
        for ii = 1: Ngroup +1
            temp = reshape(out_temp(:,ii),Nnodes,N);
            out(:,ii) = prod(temp,2);
        end
        
    end

end