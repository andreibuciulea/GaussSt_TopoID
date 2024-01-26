function reg = get_reg(model,prms)
    sig_type = prms.sig_type;
    g_type = prms.g_type;
    is_MRF = strcmp(sig_type,'MRF');
    max_iters = prms.max_iters;
    if isfield(prms,'verbose');verbose = prms.verbose; 
    else; verbose = false; end
    
    options = [g_type '-' sig_type '-' model];
    
    if strcmp(options,'ER-Poly-GSR') || strcmp(options,'ER-MRF-GSR') || strcmp(options,'ER-SSEM-GSR')
        reg = struct('epsilon',1e-6,'lambda',1e8,'rho',1,'max_iters',max_iters,'verbose',verbose);

    elseif strcmp(options,'ER-MRF-GGSR-j') || strcmp(options,'ER-MRF-GGSR-j-v1')
        reg = struct('rho',1,'delta',1e-4,'beta',1e-5,'eta',1e-4,'L1',0.01,'L2',0.2,'max_iters',max_iters);
        %reg = struct('rho',1e-2,'delta',1e-5,'beta',1e-4,'eta',1e-4,'L1',0.1,'L2',1,'max_iters',max_iters);
   
    elseif strcmp(options,'ER-SSEM-GGSR-j-v1') || strcmp(options,'BA-SSEM-GGSR-j-v1') || strcmp(options,'SBM-SSEM-GGSR-j-v1') || strcmp(options,'SW-SSEM-GGSR-j-v1')
        reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-2,'eta',1e-4,'L1',0.01,'L2',1,'max_iters',max_iters,'mit',1000);
        %reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-3,'eta',1e-4,'L1',0.01,'L2',0.1,'max_iters',max_iters,'mit',100);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(options,'ER-Poly-GGSR-j-v1') 
        reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-4,'eta',1e-4,'L1',0.01,'L2',1,'max_iters',max_iters,'mit',1000);
    elseif strcmp(options,'BA-Poly-GGSR-j-v1') 
        reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-2,'eta',1e-4,'L1',0.01,'L2',1,'max_iters',max_iters,'mit',1000);
    elseif strcmp(options,'SBM-Poly-GGSR-j-v1')
        reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-2,'eta',1e-4,'L1',0.01,'L2',1,'max_iters',max_iters,'mit',1000);
    elseif strcmp(options,'SW-Poly-GGSR-j-v1')
        reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-2,'eta',1e-4,'L1',0.01,'L2',1,'max_iters',max_iters,'mit',1000);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(options,'ER-Poly-GGSR-j')
        reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-3,'eta',1e-4,'L1',1e-1,'L2',1,'max_iters',max_iters,'mit',1000);
        %reg = struct('rho',1e-1,'delta',1e-5,'beta',1e-3,'eta',1e-3,'L1',10,'L2',10,'max_iters',max_iters);
    elseif strcmp(options,'ER-SSEM-GGSR-j')
        reg = struct('rho',1e-1,'delta',1e-4,'beta',1e-3,'eta',1e-4,'L1',1e-1,'L2',1,'max_iters',max_iters,'mit',1000);
        %reg = struct('rho',1e-1,'delta',1e-5,'beta',1e-3,'eta',1e-3,'L1',10,'L2',10,'max_iters',max_iters);

    elseif strcmp(options,'ER-MRF-GGSR-bi') 
        reg = struct('rho',1e2,'max_iters',max_iters);

    elseif strcmp(options,'ER-Poly-GGSR-bi')
        reg = struct('rho',1e2,'max_iters',max_iters);

    elseif strcmp(options,'ER-MRF-est-S-G')
        reg = struct('rho',1e-3,'lambda',1,'max_iters',max_iters,'verbose',verbose);
    elseif strcmp(options,'ER-Poly-est-S-G')
        reg = struct('rho',1e3,'lambda',1e-2,'max_iters',max_iters,'verbose',verbose);
        
    elseif strcmp(options,'ER-MRF-GGSR') || strcmp(options,'ER-SSEM-GGSR')
        reg = struct('rho',1e-2,'lambda',1e1,'tol',1e-4,'max_iters',max_iters,'verbose',verbose); %1e4-1e6   
        %reg = struct('rho',1.6e-2,'lambda',1,'max_iters',max_iters,'verbose',verbose); %1e3
        %reg = struct('rho',1e-4,'lambda',1,'max_iters',max_iters,'verbose',verbose); %1e3
    elseif strcmp(options,'ER-SSEM-GGSR') || strcmp(options,'BA-SSEM-GGSR') || strcmp(options,'SW-SSEM-GGSR') || strcmp(options,'SBM-SSEM-GGSR')
        reg = struct('rho',1e-2,'lambda',1e1,'tol',1e-4,'max_iters',max_iters,'verbose',verbose);   
    elseif strcmp(options,'ER-Poly-GGSR') || strcmp(options,'BA-Poly-GGSR') || strcmp(options,'SW-Poly-GGSR') || strcmp(options,'SBM-Poly-GGSR')     
        reg = struct('rho',1e-4,'lambda',1e2,'tol',1e-4,'max_iters',max_iters,'verbose',verbose);%sin ruido
        %reg = struct('rho',1e-2, 'lambda',1e2, 'tol',1e-4, 'max_iters',max_iters,'verbose',verbose); %con ruido
        %reg = struct('rho',1e-6,'lambda',1e3,'max_iters',max_iters,'verbose',verbose);%fig3
        %para M = 1e6, N=50 y Rho=1e-4 --> Lambda=1e4;
        %para M = 1e6, N=40 y Rho=1e-4 --> Lambda=1e4;
    elseif strcmp(options,'RBF-Poly-est-S-GST')
        reg = struct('rho',1e-5,'lambda',1e2,'max_iters',max_iters,'verbose',verbose); %M = 1e6
    elseif strcmp(options,'BA-Poly-est-S-GST')
        reg = struct('rho',1e-5,'lambda',1e2,'max_iters',max_iters,'verbose',verbose); %M = 1e6
    
    elseif strcmp(options,'ER-Poly-est-S-GST-block')
        reg = struct('rho',5e-5,'lambda',1e3,'max_iters',max_iters,'verbose',verbose); 
    elseif strcmp(options,'ER-MRF-est-S-GST-block') 
        reg = struct('rho',1e-4,'lambda',1e2,'max_iters',max_iters,'verbose',verbose);
        
    elseif strcmp(options,'ER-Poly-est-S-GST-3B') || strcmp(options,'RBF-Poly-est-S-GST-3B') || strcmp(options,'BA-Poly-est-S-GST-3B')
        reg = struct('rho',1e-5,'la1',1e4,'la2',1e3,'kappa', 1e-2,'tol',5e-2,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type); 
        %lambda = 1e4 pero para el exp5 lambda=1e5;
    elseif strcmp(options,'ER-MRF-est-S-GST-3B') || strcmp(options,'RBF-MRF-est-S-GST-3B') || strcmp(options,'BA-MRF-est-S-GST-3B')
        reg = struct('rho',1e-5,'la1',1e4,'la2',1e3,'kappa', 1e-2,'tol',1e-2,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type);
%reg = struct('rho',1e-4,'lambda',1e6,'gamma', 1,'max_iters',max_iters,'verbose',verbose);

    elseif strcmp(options,'ER-Poly-est-S-GST-3BI') || strcmp(options,'RBF-Poly-est-S-GST-3BI') || strcmp(options,'BA-Poly-est-S-GST-3BI')
        reg = struct('rho',1e-6,'lambda',1e5,'gamma',1e6,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type); 
    elseif strcmp(options,'ER-MRF-est-S-GST-3BI') || strcmp(options,'RBF-MRF-est-S-GST-3BI') || strcmp(options,'BA-MRF-est-S-GST-3BI')
        reg = struct('rho',1e-6,'lambda',1e6,'gamma', 1e6,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type);

    elseif strcmp(options,'ER-Poly-est-S-GST-3BRw') || strcmp(options,'RBF-Poly-est-S-GST-3BRw') || strcmp(options,'BA-Poly-est-S-GST-3BRw')
        reg = struct('rho',1e-6,'lambda',1e-6,'gamma',1e6,'kappa',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF); 
    elseif strcmp(options,'ER-MRF-est-S-GST-3BRw') || strcmp(options,'RBF-MRF-est-S-GST-3BRw') || strcmp(options,'BA-MRF-est-S-GST-3BRw')
        reg = struct('rho',1e-6,'lambda',1e-6,'gamma', 1e6,'kappa',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);

    elseif strcmp(options,'ER-Poly-est-S-GST-fast') || strcmp(options,'RBF-Poly-est-S-GST-fast') || strcmp(options,'BA-Poly-est-S-GST-fast')
        reg = struct('rho',8,'lambda',1e-4,'gamma',1,'kappa',1e-2,'alpha',1e-5,'tol',5e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF); 
    elseif strcmp(options,'ER-MRF-est-S-GST-fast') || strcmp(options,'RBF-MRF-est-S-GST-fast') || strcmp(options,'BA-MRF-est-S-GST-fast')
        reg = struct('rho',8,'lambda',1e-4,'gamma', 1,'kappa',1e-2,'alpha',1e-5,'tol',5e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    elseif strcmp(options,'ER-Var-est-S-GST-fast') || strcmp(options,'RBF-Var-est-S-GST-fast') || strcmp(options,'BA-MRF-est-S-GST-fast')
        reg = struct('rho',2,'lambda',1e-4,'gamma', 1,'kappa',1e-2,'alpha',1e-4,'tol',1e-4,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    
    elseif strcmp(options,'ER-Var-GGSR-fast')
        reg = struct('rho',3,'la1',1e-6,'la2',1e-6,'kappa',1e-2,'alpha',1e-4,'tol',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    elseif strcmp(options,'ER-Poly-GGSR-fast')
        reg = struct('rho',8,'la1',1e-5,'la2',1e-3,'kappa',1e-2,'alpha',1e-4,'tol',5e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    elseif strcmp(options,'ER-MRF-GGSR-fast')
        reg = struct('rho',1,'la1',1e-7,'la2',1e-3,'kappa',1e-2,'alpha',1e-4,'tol',5e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
                                                                                         %5e-2   
    elseif strcmp(options,'ER-Poly-est-S-GST-3BFast') || strcmp(options,'RBF-Poly-est-S-GST-3BFast') || strcmp(options,'BA-Poly-est-S-GST-3BFast')
        reg = struct('rho',1e-6,'lambda',1e5,'gamma',1e6,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF); 
    elseif strcmp(options,'ER-MRF-est-S-GST-3BFast') || strcmp(options,'RBF-MRF-est-S-GST-3BFast') || strcmp(options,'BA-MRF-est-S-GST-3BFast')
        reg = struct('rho',1e-6,'lambda',1e6,'gamma', 1e6,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    
    elseif strcmp(options,'ER-MRF-GL') || strcmp(options,'ER-Poly-GL') 
        reg = struct('rho',1e-6,'tol',1e4,'max_iters',max_iters*100,'verbose',verbose);
    elseif strcmp(options,'ER-SSEM-GL')
        reg = struct('rho',1e-6,'tol',1e4,'max_iters',max_iters*100,'verbose',verbose);

    %real data regularizers for financial
    elseif strcmp(model,'GSR-finance')
        reg = struct('epsilon',1e-6,'lambda',1e8,'max_iters',max_iters,'verbose',verbose);
    elseif strcmp(model,'est-S-G-finance')
        reg = struct('rho',1e-4,'lambda',1e3,'max_iters',max_iters,'verbose',verbose);
    elseif strcmp(model,'est-S-GST-finance')
        reg = struct('rho',1e3,'lambda',1e-3,'max_iters',max_iters,'verbose',verbose);
        sig_type = 'MRF';
    elseif strcmp(model,'est-S-GST-sparse-finance')
        reg = struct('rho',1e-2,'lambda',1.2e-1,'max_iters',max_iters,'verbose',verbose);
        sig_type = 'MRF';
        
    elseif strcmp(model,'GL-finance')
        reg = struct('rho',0.1,'tol',1e5,'max_iters',1e5,'verbose',verbose);
    
    elseif strcmp(model,'est-S-GST-3B-finance')
        reg = struct('rho',1e2,'lambda',1e-3,'gamma', 1e1,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    elseif strcmp(model,'est-S-GST-3BRw-finance')
        reg = struct('rho',1e-6,'lambda',1e-6,'gamma', 1e6,'kappa',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    elseif strcmp(model,'est-S-GST-fast-finance')
       reg = struct('rho',8,'lambda',1e-4,'gamma',1,'kappa',1e-2,'alpha',1e-5,'tol',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    elseif strcmp(model,'GST-fast-finance')
       reg = struct('rho',1,'la1',1e-4,'la2',1e-4,'gamma',1,'kappa',1e-2,'alpha',1e-4,'tol',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
       %reg = struct('rho',1,'la1',1e-3,'la2',1e-3,'gamma',1,'kappa',1e-2,'alpha',1e-4,'tol',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);

    %real data regularizers for proteins
    elseif strcmp(model,'GSR-noise-prot')
        reg = struct('epsilon',1.15,'max_iters',max_iters,'verbose',verbose);
    elseif strcmp(model,'est-S-G-prot')
        reg = struct('rho',1e-4,'lambda',1e3,'max_iters',max_iters,'verbose',verbose);
    elseif strcmp(model,'est-S-GST-prot')
        reg = struct('rho',0.42,'lambda',770,'max_iters',max_iters,'verbose',verbose);
    elseif strcmp(model,'est-S-GST-noise-prot')
        reg = struct('rho',1e-5,'lambda',1e2,'epsilon',1,'max_iters',max_iters,'verbose',verbose);
        reg.isMRF = true;
    elseif strcmp(model,'GL-prot')
        reg = struct('rho',1e-6,'tol',1e4,'max_iters',max_iters*100,'verbose',verbose);
    elseif strcmp(model,'est-S-GST-3B-noise-prot')
        %reg = struct('rho',5.3e-3,'lambda',5.6234e3,'max_iters',max_iters,'verbose',verbose); 
        reg = struct('rho',7.4e3,'lambda',1e3,'max_iters',max_iters,'verbose',verbose);
    elseif strcmp(model,'est-S-GST-3B-prot')
        reg = struct('rho',1e2,'lambda',0.1,'gamma', 1e2,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type);
    elseif strcmp(model,'est-S-GST-3BRw-prot')
        %reg = struct('rho',0.08,'lambda',0.1,'gamma', 1e2,'kappa',1e-2,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type);
        reg = struct('rho',1e-3,'lambda',1e-1,'gamma', 1e3,'kappa',1e-4,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type);

    %regularizers for real data senate
    elseif strcmp(model,'GSR-senate')
        reg = struct('lambda',1e3,'max_iters',max_iters,'verbose',verbose);
    
    elseif strcmp(model,'GL-senate')
        reg = struct('rho',0.8,'tol',1e4,'max_iters',max_iters*100,'verbose',verbose);
    
    elseif strcmp(model,'est-S-GST-3BRw-senate')
        reg = struct('rho',1e-2,'lambda',1e-2,'gamma', 1e4,'kappa',1e-2,'max_iters',max_iters,'verbose',verbose,'C_type',sig_type,'is_MRF',is_MRF);
    elseif strcmp(model,'est-S-GST-fast-senate')
       reg = struct('rho',1.4,'lambda',1e-6,'alpha',1e-6,'gamma',1,'kappa',1e-2,'tol',1e-2,'max_iters',max_iters,'verbose',verbose,'is_MRF',is_MRF);
    else
        reg = {};
        error(['Error: unknown reg for ' options])
    end

    reg.isMRF = strcmp(sig_type,'MRF'); %cuidado para datos reales
    
        
end