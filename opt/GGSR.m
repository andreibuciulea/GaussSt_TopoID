function [S_hat,out] = GGSR(C,regs)
    C = C/sqrt(max(eig(C)));

    N = size(C,1);
    rho = regs.rho;
    lambda = regs.lambda;
    max_iters = regs.max_iters;
    %max_iters = 20;
    verbose = regs.verbose;
    tol = regs.tol;
    %isMRF = regs.is_MRF;
    S_prev = zeros(N);%C-diag(diag(C));
    Pr_prev = eye(N);%inv(C);
    S_hat = S_prev;
    O_F = zeros(max_iters+1,1);

    errS = zeros(max_iters,1);
    e_time = zeros(max_iters,1);
    fobj = zeros(max_iters,1);
    S_true = regs.S_true;
    nSt = norm(S_true,'fro')^2;
    my_time = tic;

    %lambda = 10;
    %rho = 1e-6;
    %figure()
    for i = 1:max_iters
        %Subproblem I 
        cvx_begin quiet
            variable Pr(N,N) symmetric semidefinite

            minimize(-log_det(Pr)+trace(C*Pr)+lambda*norm(Pr*S_hat-S_hat*Pr,'fro'))
            
            subject to
                %if isMRF
                %    diag(Pr) >= 0;
                %    Pr*ones(N,1) >= 1;
                %end
        cvx_end
        
        if sum(sum(isnan(Pr))) == 0
           Pr_prev = Pr;
        else
            Pr = Pr_prev;
        end
        if verbose
            disp(['Problem 1: logdet:' num2str(-log_det(Pr)) '  trace(C*Pr):' num2str(trace(C*Pr))  '  Fro:' num2str(norm(Pr*S_hat-S_hat*Pr,'fro'))])
        end

        %Subproblem II

        cvx_begin quiet
            variable S_hat(N,N) symmetric

            minimize(rho*norm(S_hat(:),1)+lambda*norm(Pr*S_hat-S_hat*Pr,'fro'))

            subject to
                diag(S_hat) <= 1e-6;
                S_hat >= 0;
                S_hat*ones(N,1) >= 1;
        cvx_end
        %S_hat = S_hat/max(max(S_hat));
        if verbose
            disp(['Problem 2: norm(S,1):' num2str(norm(S_hat(:),1)) '  Fro Comm:' num2str(norm(Pr*S_hat-S_hat*Pr,'fro')) '  Fro S:' num2str(norm(S_hat-S_prev,'fro')^2)])
        end
        %stop criteria
        fobj(i) = obj_fun(C,Pr,S_hat,rho,lambda);
        obj_f = objective_function(Pr, S_hat, C, rho, lambda, zeros(N));
        O_F(i+1) = obj_f.value;
        errS(i) = norm(S_hat-S_true,'fro')^2/nSt;
        e_time(i) = toc(my_time);

        if abs(O_F(i)-O_F(i+1)) <= tol
            break
        end
        
        if sum(sum(isnan(S_hat))) == 0
           S_prev = S_hat;
        else
            S_hat = S_prev;
        end


        if true
            figure(1)
            subplot(121)
            imagesc(Pr)
            title('Pr')
            colorbar()
            subplot(122)
            imagesc(S_hat)
            title('S hat')
            colorbar()   
        end

    end
    out.objective = O_F;
    out.S_hat = S_hat;
    out.Pr = Pr;
    out.errS = errS;
    out.e_time = e_time;
    out.fobj = fobj;
end