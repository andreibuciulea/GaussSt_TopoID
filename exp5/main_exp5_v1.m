clear all
addpath('../utils');
addpath('../opt');
addpath(genpath('../../global_utils'));

%rng(4)

g_type = 'ER'; 
N = 20;          
M = 1e6;
sig_type = 'Poly';     Ct = numel(sig_type);
Models = {'GL','GSR','GGSR','GGSR-jiaxi-v1'}; %models = {'gti-GST'}; 
Models = {'GGSR','GGSR-jiaxi'}; %models = {'gti-GST'}; 
nM = numel(Models);
sG = [20,30,40,50,60,70,80];
%sG = [20,30,40,50,60,70,80,90,100,110,120,130,140,150];
nsG = numel(sG);
nG = 64;
mit = 1000;

%%%Graph gen params
    %parameters for graph generation
    p = 0.1;
%%%Signal gen params 
    max_iters = 10;
    sampled = true;
    L = 3;

fsc_S = zeros(nG,nsG);
err_T = zeros(nG,nsG);
err_S = zeros(nG,nsG);
est_obj_fun = zeros(nG,nsG);
true_obj_fun = zeros(nG,nsG);

errS = cell(nG,1);%zeros(nG,nI,max_iters);
fobj = cell(nG,1);%zeros(nG,nI,max_iters);
e_time = cell(nG,1);%zeros(nG,nI,max_iters);
t0 = tic;
parfor ng = 1:nG %number of graphs
    errSi = cell(nsG,nM);%zeros(nI,max_iters);
    fobji = cell(nsG,nM);%zeros(nI,max_iters);
    e_timei = cell(nsG,nM);%= zeros(nI,max_iters);

    disp(['Graph: ' num2str(ng)])
    
    %general parameters
    prms = struct('sampled',sampled,'max_iters',max_iters,...
                  'M',M,'L',L,'g_type',g_type,'sig_type',sig_type);
    for ns = 1:nsG
        %toc(t0)
        %disp(['Graph size: ' num2str(sG(ns))])
        N = sG(ns);
        g_prms = struct('N',N,'ER_p',p,'g_type',g_type);   
    
        %Generate graphs
        S = generate_graph(g_prms).A;
        %while min(abs(eig(eye(N)-S))) <= 1e-1  
            %disp(min(abs(eig(eye(N)-S))))
       %     S = generate_graph(g_prms).A;
       % end
    
        %Generate signals
        out = generate_graph_signals(S, prms);
        C = out.C;
        Theta = out.C_inv;
        for m = 1:nM
            %get regularizers and add S and Theta true
            model = Models{m};
            regs = get_reg(model,prms);
            regs.S_true = S;
            regs.Theta = Theta;
            regs.mit = mit; 
    
            true_obj_fun(ng,ns,m) = obj_fun(C,Theta,S,regs.rho,1e-4);
           
            [S_hat,out] = estimate_S(C,model,regs);
            S_hat = out.S_hat;
            %remove the iteration limit
            %The value of error in S for each iteration
            errSi{ns,m} = out.errS;
            %The value of the Objective Function for each iteration
            fobji{ns,m} = out.fobj;
            %The elapsed time for each iteration
            e_timei{ns,m} = out.e_time;
    
            S_hat = S_hat/max(max(S_hat));
            est_obj_fun(ng,ns,m) = obj_fun(C,out.Pr,S_hat,regs.rho,1e-4);
            Pr = out.Pr/max(max(out.Pr));
            %fsc_S(ng,ni) = fscore(S,S_hat);
            err_S(ng,ns,m) = norm(S - S_hat,'fro')^2/norm(S,'fro')^2;
            err_T(ng,ns,m) = norm(Theta-out.Pr,'fro')^2/norm(Theta,'fro')^2;
        end
    end
    %The value of error in S for each iteration
    errS{ng} = errSi;
    %The value of the Objective Function for each iteration
    fobj{ng} = fobji;
    %The elapsed time for each iteration
    e_time{ng} = e_timei;
end 
toc(t0)

%%
% figure()
% plot(true_obj_fun)
% grid on
% hold on
% plot(est_obj_fun)
% grid on
% legend('True','Est')
mrkt = {'-',':';'-',':';'-',':';'-',':';'-',':';'-',':';'-',':';'-',':';'-',':';'-',':';'-',':';'-',':'};
lgd = {'20-GGSR','20-GGSR-J','30-GGSR','30-GGSR-J',...
    '40-GGSR','40-GGSR-J','50-GGSR','50-GGSR-J'};
%lgd = {'20-GGSR','20-GGSR-J','30-GGSR','30-GGSR-J','40-GGSR','40-GGSR-J'};
for g = 1:10
    figure()
    for ni = 1:nsG
      for nm = 1:nM   
         subplot(211)
%         semilogy(errS{g}{ni})
%         title('Error in S')
%         %legend(text_nI{1:nI})
%         hold on
%         grid on
%         xlabel('Number of iterations')
%         subplot(223)
%         semilogy(fobj{g}{ni})
%         title('Objective function')
%         %legend(text_nI{1:nI})
%         grid on
%         hold on
%         xlabel('Number of iterations')
        %subplot(222)
        loglog(e_time{g}{ni,nm}, errS{g}{ni,nm},mrkt{ni,nm},'Markersize',0.1,LineWidth=2)
        hold on
        %legend(lgd)
        title('Error in S')
        grid on
        xlabel('Time (s)')

        subplot(212)
        loglog(e_time{g}{ni,nm}, fobj{g}{ni,nm},mrkt{ni,nm},'Markersize',0.1,LineWidth=2)
        hold on
        %legend(lgd)
        title('Objective function')
        grid on
        xlabel('Time (s)')
%         loglog(e_time{g}{ni,nm},fobj{g}{ni,nm})
%         hold on
%         title('Objective function')
%         %legend(text_nI{1:nI})
%         grid on
%         xlabel('Time (s)')
      end
    end

end


%%
all_time = zeros(nG,nsG,nM);
all_err = zeros(nG,nsG,nM);
for g = 1:nG
    for ni = 1:nsG
      for nm = 1:nM  
        all_time(g,ni,nm) = max(e_time{g}{ni,nm});
        aux = errS{g}{ni,nm};
        all_err(g,ni,nm) = aux(end);
      end
    end
end

figure()
semilogy(sG,squeeze(median(all_time)))
grid on

figure()
semilogy(sG,squeeze(median(err_S)))
grid on

%%
all_time = zeros(nG,nsG,nM);
all_err = zeros(nG,nsG,nM);
for g = 1:nG
    for ni = 1:nsG
      for nm = 1:nM  
        aux = errS{g}{ni,nm};
        all_errs(g,ni,nm) = aux(end);
      end
    end
end

figure()
semilogy(sG,squeeze(median(all_time)))
grid on
