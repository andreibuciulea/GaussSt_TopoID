function [S_out,out] = GGSR_fast(C,regs)

    %C = C/max(max(C));
    
    N = size(C,1);
    la1=1e-2;%regs.la1; %SPR I stationarity 
    la2=1e-3;%regs.la2;
    la3=1e-3;
    
    rho = regs.rho;
    alpha = regs.alpha;
    kappa = regs.kappa;
    tol = regs.tol;
    max_iters = regs.max_iters;
    %S_true = mbinarize(C-diag(diag(C)),2);
    S_true = regs.S_true;
    aux = find(S_true(:,1)~=0);
    n2 = aux(1);

    verbose = regs.verbose;
    
    kI = kappa*eye(N);
    %C = C + kI;
    
    
    Pr2 = eye(N);
    %Pr2 = diag(1./diag(C)); 
    Z = zeros(N);
    Y = zeros(N);
    S = zeros(N);
    %S=regs.S_true;
    a = 1;
    mit = 1e3;
    S_prev = zeros(N);
    spr = nan(3,max_iters);
    O_F = zeros(max_iters+1,1);
    for t = 1:max_iters
        t1 = tic;
        %%%%%%%%%%%%%SubproblemII%%%%%%%%%%%%%%%
        %estimate Pr1 given Pr2, S, Y, Z. 
        [V,D] = eig(S);
        d = diag(D);
        E1 = la2*Pr2-C+S*Y'-Y'*S+Z';
        E = V'*E1*V;
        D_aux = repmat(d.^2,1,N);
        X_aux = E./(la2*ones(N)+la3*D_aux+la3*D_aux'-2*la3*(d*d'));
        Pr1 = V*X_aux*V';
        
        A = Pr1*S-S*Pr1;
        B = Pr2-Pr1;
        spr(1,t)=toc(t1);

        if verbose 
            disp(['Problem1: Tr(C*Pr1): ' num2str(trace(C*Pr1)) '  ||A||fro: ' num2str(norm(A,'fro')) '  ||B||fro: ' num2str(norm(B,'fro'))...
            '  Tr(Y*A): ' num2str(trace(Y*A)) '  Tr(Z*B): ' num2str(trace(Z*Pr1))])
        end

        t2 = tic;
        
        %%%%%%%%%%%%%SubproblemII%%%%%%%%%%%%%%%
        %estimate Pr2 given Pr1, and Z. 
        [U,Gamm] = eig(la2/a*(Pr1+kI)-Z/a);
        Pr2 = 1/(2*la2/a)*U*(Gamm+sqrt(Gamm.^2 + 4*la2/a*eye(N)))*U'-kI;
        %Pr2 = 1/la2*U*(Gamm+sqrt(Gamm.^2 + 4/la2*eye(N)))*U'-kI;
        spr(2,t)=toc(t2);
        if verbose
            disp(['Problem2: LogDet(Pr2): ' num2str(-log(det(Pr2))) '  Tr(Z*Pr2): ' num2str(trace(Z*Pr2)) '  ||Pr1-Pr2||fro: ' num2str(norm(Pr1-Pr2,'fro'))])
        end
        %Pr2(Pr2<0.05) = 0;
        t3 = tic;
        %%%%%%%%%%%%%SubproblemIII%%%%%%%%%%%%%%%
        %estimate S given Pr2, and Y.
        la_max = max(eig(Pr2))^2;
        gamma = rho/(4*la1*la_max);
        %gamma = 1;
        for ii = 1:mit
            grad_S = la1*((S*Pr2-Pr2*S)*Pr2-Pr2*(S*Pr2-Pr2*S)) + Pr2*Y'-Y'*Pr2; %Y*Pr2-;
            Dst = S-gamma*grad_S; 
            S_hat = max(0,Dst-diag(diag(Dst))-alpha);
            S_hat(n2,2) = 1;
            S_hat(2,n2) = 1;
            S = S_hat;
            S = S/max(max(S));
        end
        
        spr(3,t) = toc(t3);


        if verbose
            disp(['Problem3: ||S||1: ' num2str(norm(S(:),1)) '  ||A||fro: ' num2str(norm(A,'fro')) '  Tr(Y*A): ' num2str(trace(Y*A))])
            fprintf('\n');
        end

        obj_f = objective_function(Pr2, S, C, rho, la1, zeros(N));
        O_F(t+1) = obj_f.value;
        disp(abs(O_F(t)-O_F(t+1)))
        if abs(O_F(t)-O_F(t+1)) <= tol
            figure()
            plot(O_F)
            break
        end

        %uptade Z and Y
        Z = Z+la2*(Pr2-Pr1);
        Y = Y+la1*(S*Pr2-Pr2*S);
       
        if verbose
            figure(2)
            set(gca,'units','normalized','outerposition',[0 0 1 1])
            subplot(231)
            imagesc(Pr1)
            title('Pr1')
            colorbar()
            subplot(232)
            imagesc(Pr2)
            title('Pr2')
            colorbar()
            subplot(233)
            imagesc(S)
            title('S')
            colorbar()
            subplot(234)
            imagesc(Z)
            title('Z')
            colorbar()
            subplot(235)
            imagesc(Y)
            title('Y')
            colorbar()
            subplot(236)
            imagesc(C)
            colorbar()
            title(['C'])
            disp(['S frob:' num2str(norm(S-S_prev,'fro'))])   
        end
        out.spr= spr;
        S_prev = S;
    end
    %Fcur = objective_function(Pr2, S, C, rho, la1);
    out.objective = O_F;
    S_out = S; 
    out.S = S;
    out.Theta = Pr2;
    out.Pr = Pr2;
    
end