function out = GGSR_j(C,prms)
    
    %This algorithm is the same as GGSR_jiaxi adding obj function and
    %estimation error at each step

    if isfield(prms,'constraint'); constraint = prms.constraint; 
        else; constraint = 0; prms.constraint = 0; end
    if isfield(prms,'max_iters'); max_iters = prms.max_iters; 
        else; max_iters = 50; prms.max_iters = 50; end
    if isfield(prms,'delta'); delta = prms.delta; 
        else; delta = 0.1; prms.delta = 0.1; end
    if isfield(prms,'rho'); rho = prms.rho; 
        else; rho = 0.1; prms.rho = 0.1; end
    if isfield(prms,'beta'); beta = prms.beta; 
        else; beta = 0.1; prms.beta = 0.1; end
    if isfield(prms,'eta'); eta = prms.eta; 
        else; eta = 0.1; prms.eta = 0.1; end
    if isfield(prms,'L1'); L1 = prms.L1; 
        else; L1 = 0.1; prms.L1 = 0.1; end
    if isfield(prms,'L2'); L2 = prms.L2; 
        else; L2 = 0.1; prms.L2 = 0.1; end
    if isfield(prms,'mit'); mit = prms.mit; 
        else; mit = 1e3; prms.mit = 1e3; end

%     if mit == 1
%         max_iters = max_iters*1000;
%     elseif mit == 10
%         max_iters = max_iters*100;
%     elseif mit == 100
%         max_iters = max_iters*10;
%     elseif mit == 1000
%         max_iters = max_iters*1;
%     end
    N = size(C,1);

    P  = zeros(N);Z = P;Q = P;Y = P;
    T = eye(N);
    S = zeros(N);
    
    fobj = zeros(max_iters,1);
    dobj = zeros(max_iters,1);
    errS = zeros(max_iters,1);
    e_time = zeros(max_iters,1);


    S_true = prms.S_true;
    nSt = norm(S_true,'fro')^2;
    my_time = tic;


%     prms.beta = 1e-2;
%     prms.L1 = 1;
%     prms.L2 = 0.1;
%     prms.rho = 1;
%     prms.mit = 100;

    for t = 1:max_iters
        if beta > 1e3
            beta = beta*(beta/(beta+0.05*beta));
            prms.beta = beta;
        end
    %%%%% For the second subproblem we update T, Q, and Y 
        [T,Q,Y] = gti_T(C,S,T,Q,Y,prms);

    %%%%% For the first subproblem we update S, P, and Z
        [S,Z,P] = gti_S(T,S,Z,P,prms);
        
        %S = S/max(max(S));
        fobj(t) = obj_fun(C,T,S,prms.rho,prms.beta);
        errS(t) = norm(S-S_true,'fro')^2/nSt;
        e_time(t) = toc(my_time);
        if t > 1
            dobj(t) = abs(fobj(t-1)-fobj(t)); 
            if dobj(t) < 5e-3
               %disp(['Converged t=' num2str(t) ' time: ' num2str(e_time(t))])
               %break; 
            end
        end
        if mod(t,1) == 0
            figure(3)
            subplot(121)
            imagesc(S)
            colorbar()
            subplot(122)
            imagesc(T)
            colorbar()
            figure(2)
            subplot(211)
            semilogy(errS)
            grid on
            title('Value of the estiamtion error for S')
            xlim([0 max_iters])
            subplot(212)
            semilogy(dobj)
            grid on
            title('Value of the bjective function')
            xlim([0 max_iters])
        end
    end

    out.Pr = T;
    out.S_hat = S;
    out.fobj = fobj;
    out.errS = errS;
    out.e_time = e_time;
end