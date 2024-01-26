function [S_hat, out]= GSR(C,reg)

N = size(C,1);

verbose = reg.verbose;
lambda = reg.lambda;
max_iters = 50;%reg.max_iters;
epsilon = 1e-3;%reg.epsilon;
%max_iters = 15;

if verbose
   disp('  -Starting GSR Low Rank optimization...') 
end

for i = 1:max_iters
    cvx_begin quiet
        variable S_hat(N,N) symmetric
        minimize (norm(S_hat(:),1)) %+ lambda/2*square_pos(norm(C*S_hat - S_hat*C, 'fro')))
        subject to
            abs(diag(S_hat)) <= 1e-6;
            S_hat >= 0;
            S_hat*ones(N,1) >= 1;
            norm(C*S_hat - S_hat*C, 'fro') <= epsilon;
    cvx_end
    if strcmp(cvx_status,'Solved') %%|| strcmp(cvx_status,'Failed')  
        if verbose
            disp(['Break after ' num2str(i) ' iterations'])
        end
        break;
    else
        if verbose
            disp('Updating epsilon ...')
        end
        epsilon = epsilon*1.5;
        %epsilon = epsilon*2;
    end
end
out.S_hat = S_hat;
out.Pr = diag(sum(S_hat))-S_hat;
out.errS = 0;
out.fobj = 0;
out.e_time = 0;
end