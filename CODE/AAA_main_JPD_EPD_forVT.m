
% This code illustrates how to obtain the systemic risk measures such as 
% JPD and EPD given the estimation of a skew t-t factor copula with GAS 
% dynamics on a sample of 10 variables. Although this code assumes that 
%  all estimation is done by "variacbe targetting" method based on 
% "AAA_main_GASFacCop_Skewtt_VT.m", it is easily modified to accommodate the
% estimation by non-VT method based on "AAA_main_GASFacCop_Skewtt_Ngroup.m"
%  
%
%  Dong Hwan Oh and Andrew Patton
%
%  16 February 2015
%
%  This code is to accompany the paper:
%  Oh, D.H. and A.J. Patton, 2015, Time-Varying Systemic Risk: Evidence 
%  from a Dynamic Copula Model of CDS Spreads, working paper, Duke University. 

input_path = 'D:\Core\Matlab_programs\Factor_copulas\GAS_factor_copula_toolbox\';  % where the copula data was saved
output_path = 'D:\Core\Matlab_programs\Factor_copulas\GAS_factor_copula_toolbox\'; % where the results MAT file should be saved

load([input_path,'sample_CDS_data.mat'],'sample_CDS_data');   % CDS spreads data
load([input_path,'sample_data.mat'],'copula_data');           % U's obtained from the above CDS spreads 

load([input_path,'results_VT_from_GASFactorToolbox2.mat'],'theta_NMLE','lam_bar','lam');    % estimation results from VT 
load([input_path,'sample_marignal_data.mat'],'cond_mean','cond_var','eachfirm_resid','eachfirm_return','GARCH_para_all','AR_para_all') %estimation results from marginal dist.


Nsim = 1000 ;   % # of replications 
T_1year = 250;  % # of days for one year

[TT, NN] = size(copula_data);

group_code = [1 2 3 4 5 6 7 8 9 10] ;
Ngroup = max(group_code);

monthly_date = (1:20:TT);   %50 months

alpha     = theta_NMLE(1);
beta      = theta_NMLE(2);
nuinv_z   = theta_NMLE(3);
nuinv_eps = theta_NMLE(4);
psi_z     = theta_NMLE(5);

omega = (1-beta) * log(lam_bar) ;

%%%
%%% Nodes and weight for Gauss-Legandre quadrature 
Npt_Quad = 50;              
[Wa, Wb] = GLNodeWt(Npt_Quad);
GLweight = [Wa, Wb];


%%% Evaluate the marginal cdf (G) and pdf (g) of skew t - t factor model at various x and factor loadings fixing other
% parameters such as nuinv_z, nuinv_eps, and psi_z.
x1 = -15;
x2 = 15;
Npoints = 100;
x_grid = (x1:(x2-x1)/(Npoints-1):x2)';
x_grid = [-30;x_grid;30];

lam_grid = [0.001,0.01,(0.05:0.05:2.5)]' ;

[Gcdf_ini, Gpdf_ini] = Gcdfpdf_Skewtt([nuinv_z; nuinv_eps; psi_z], GLweight, x_grid, lam_grid);

%%% Interpolate the marginal cdf and pdf along with finer grids of factor loadings (lam)
dense_lam_grid = [0.001;(0.01:0.001:2.5)'];
Gcdf = nan(length(x_grid),length(dense_lam_grid));
Gpdf = nan(length(x_grid),length(dense_lam_grid));

for i = 1:length(dense_lam_grid)
    Gcdf(:,i) = interp2(lam_grid, x_grid, Gcdf_ini, dense_lam_grid(i), x_grid);
    Gpdf(:,i) = interp2(lam_grid, x_grid, Gpdf_ini, dense_lam_grid(i), x_grid);
end


%%% Simulate 1-year-ahead (250-day-ahead) CDS spreads 'Nsim' times

%%% The simulation below takes long time. We recommend to break the for-loop below 
%%% (for tt_montly = ...) into several ones and utilize some parallel computing. 

Sdata_250dayNSt = nan(T_1year, NN, Nsim, length(monthly_date));

for tt_monthly =  monthly_date 
    
    eachfirm_return_ss = nan(T_1year, NN, Nsim);
    
    for ss = 1 : Nsim
        
        %%%%%%%%%%%% initializing variables %%%%%%%%%%%%%%%%
        sim_data_beforeCDF = nan(T_1year,NN);
        U    = nan(T_1year,NN);
        eta  = nan(T_1year,NN);
        
        s     = nan(T_1year,Ngroup);
        
        lam_1year = nan(T_1year+1,Ngroup);
        log_lam_1year = nan(T_1year+1,Ngroup);
                
        eachfirm_resid_ss  = nan(T_1year,NN);
        cond_var_ss        = nan(T_1year,NN);
        cond_mean_ss       = nan(T_1year,NN);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rand('state',123456*(tt_monthly) + ss  );
        u1 = rand(T_1year,1);
        
        common_factor = skewtdis_inv(u1,1/(nuinv_z),psi_z);
        idio_syc = nan(T_1year,NN);
        
        for nn = 1: NN
            rand('state',123456*(tt_monthly) + 4321*ss + nn  );
            u2 = rand(T_1year,1);
            idio_syc(:,nn) = skewtdis_inv(u2,1/(nuinv_eps),0);
        end
        
        lam_1year(1,:) = lam(tt_monthly,:);
        log_lam_1year(1,:) = log(lam_1year(1,:));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% simulate copula U's at each t and GAS updating  %%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for tt_1yr = 1: T_1year
            
            epsi = 0.001 ;
            
            if tt_1yr ~= 1
                log_lam_1year(tt_1yr,:) = omega' + alpha*s(tt_1yr-1,:) + beta*log_lam_1year(tt_1yr-1,:);
                lam_1year(tt_1yr,:) = exp(log_lam_1year(tt_1yr,:));
            end
            
            lam_1year(tt_1yr, (lam_1year(tt_1yr,:) < 0.01)) = 0.01;
            lam_1year(tt_1yr, (lam_1year(tt_1yr,:) > 2.5 )) = 2.5;
            
            for gg = 1: Ngroup
                [~, inx_a] = min( abs( lam_grid - lam_1year(tt_1yr,gg) ) );
                G_appro = Gcdf(:,inx_a);
                inx_group = find(group_code == gg);
                
                sim_data_beforeCDF(tt_1yr, inx_group) = lam_1year(tt_1yr,gg)*common_factor(tt_1yr) + idio_syc(tt_1yr,inx_group);
                U(tt_1yr,inx_group)  = spline(x_grid,G_appro,sim_data_beforeCDF(tt_1yr,inx_group)');
            end
            
            [~, N_derivative] = LL_eacht_GASFacCop_Skewtt_Ngroup([lam_1year(tt_1yr,:)'; nuinv_z; nuinv_eps; psi_z], U(tt_1yr,:)', GLweight, Gcdf, Gpdf, x_grid, dense_lam_grid, group_code, epsi) ;
            s(tt_1yr,:) = N_derivative' .* lam_1year(tt_1yr,:) ;
            
        end
        
        %%%%%%  skew t-inv (eta) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for nn = 1:NN
            
            nu_skewt  = 5 ;     % replace this line with a estimated parameter of each marginal dist.
            psi_skewt = 0.5 ;   % replace this line with a estimated parameter of each marginal dist.
            
            eta(:,nn)  = skewtdis_inv(U(:,nn) , nu_skewt, psi_skewt) ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%% con_mean  & con_var %%%%%%%%%%%%%%%%%%%%%%%%
        
        eachfirm_return_ss(1,:,ss)   = eachfirm_return(tt_monthly,:);
        eachfirm_resid_ss(1,:)       = eachfirm_resid(tt_monthly,:);
        cond_mean_ss(1,:) = cond_mean(tt_monthly,:) ;
        cond_var_ss(1,:)  = cond_var(tt_monthly,:)  ;
        
        
        for tt = 2: T_1year
            
            for nn = 1: NN
                AR_para = AR_para_all(:,nn);
                cond_mean_ss(tt,:) = AR_para'*[1; eachfirm_return_ss(tt-1,nn,ss) ];
                
                GARCH_para = GARCH_para_all(:,nn);
                cond_var_ss(tt,:)  = GARCH_para'*[1; eachfirm_resid_ss(tt-1,nn)^2 ; cond_var_ss(tt-1,nn) ];
            end
            
            eachfirm_return_ss(tt,:,ss) = cond_mean_ss(tt,:) + sqrt(cond_var_ss(tt,:)).*eta(tt,:) ;
            eachfirm_resid_ss(tt,:) = sqrt(cond_var_ss(tt,:)).*eta(tt,:);
            
        end
        
    end
    inx_temp = find(monthly_date == tt_monthly);
    Sdata_250dayNSt(:,:,:,inx_temp) = eachfirm_return_ss;   
   
end

%save('Sim_data_1yr_ahead_monthly.mat','Sdata_250dayNSt')
%load('Sim_data_1yr_ahead_monthly.mat','Sdata_250dayNSt')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% JPD & EPD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JPD = nan(length(monthly_date), 4);
EPD = nan(length(monthly_date), NN);

for tt = 1: length(monthly_date)
    
    tt_monthly = monthly_date(tt);
    
    CDS_spread_att = sample_CDS_data(tt_monthly,:);
    sim_data = Sdata_250dayNSt(:,:,:,tt);
    
    % Nsim by NN %
    Sdata_spread_250 = CDS_spread_att(ones(Nsim,1),:).*exp((squeeze(sum(sim_data,1)))') ; % 250 days ahead
    
    mean_crv99   =  3.0251 ; % replace this number with your own threshold;
    
    D_it250 = (Sdata_spread_250 >= mean_crv99) ;
    sum_across_firm = sum(D_it250, 2);
        
    Bnum = [2 3 4 5];
    for k = 1 : length(Bnum)
        
        JPD(tt,k) = mean(sum_across_firm > Bnum(k))  ;
        
    end
    
    for nn = 1:NN
        inx_row = (D_it250(:,nn) == 1);
        if sum(inx_row) == 0            
            EPD(tt,nn) =  0;
        else
            EPD(tt,nn) = mean(sum(D_it250(inx_row,:),2));
        end
    end
    
end

figure
plot(JPD)

figure
plot(EPD)

