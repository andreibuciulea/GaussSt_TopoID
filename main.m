clear all
addpath('utils');
addpath('global_utils')
addpath('opt');


rng(2)

g_type = {'ER'};        Gt = numel(g_type);
g_size = [20];          Gs = numel(g_size);
%n_samples = [1e2,1e3,1e4,1e5,1e6]; 
n_samples = [1e6];
nS = numel(n_samples);
sig_type = {'MRF'};     Ct = numel(sig_type);
models = {'GSR','GGSR','GGSR-fast','GGSR-bi','GGSR-j'}; 
models = {'GGSR','GGSR-j'}; 
nA = numel(models);

nG = 64;

%%%Graph gen params
    %parameters for graph generation
    p = 0.1;
%%%Signal gen params 
    max_iters = 50;
    sampled = true;
    L = 5;

fsc_S = zeros(nG,Gt,Gs,nS,Ct,nA);
err_T = zeros(nG,Gt,Gs,nS,Ct,nA);
err_S = zeros(nG,Gt,Gs,nS,Ct,nA);
all_spr = zeros(nG,3,max_iters);
all_differ = zeros(nG,max_iters,1e4);
t0 = tic;
parfor ng = 1:nG %number of graphs
    disp(['Graph: ' num2str(ng)])
    %initialize variables
    fsc_Sg = zeros(Gt,Gs,nS,Ct,nA);
    err_Tg = zeros(Gt,Gs,nS,Ct,nA);
    err_Sg = zeros(Gt,Gs,nS,Ct,nA);
    %general parameters
    prms = struct('sampled',sampled,'max_iters',max_iters,...
                  'L',L);
    g_prms = struct('p',p);     
    for gt = 1:Gt %number of graph types
        g_prms.g_type = g_type{gt}; 
        prms.g_type = g_type{gt};
        for gs = 1:Gs %number of graph sizes
            g_prms.N = g_size(gs);
            g_prms.p = p/(g_prms.N/20);
            S = generate_graph(g_prms).A;  %S = S/norm(S,'fro');
            for ns = 1:nS %number of samples
                prms.M = n_samples(ns);
                for ct = 1:Ct %number of types of C
                    prms.sig_type = sig_type{ct};
                    out = generate_graph_signals(S, prms);
                    C = out.C;
                    Theta = out.C_inv;
                    %Theta = Theta/max(max(Theta));
                    %C = C/max(abs(eig(C)));
                    for na = 1:nA %number of algorithms
                        model = models{na};
                        regs = get_reg(model,prms);
                        regs.S_true = S;
                        regs.Theta = Theta;
                        disp(['-logdet(T):' num2str(log(det(Theta))) ...
                              ' tr(C*T):' num2str(trace(C*Theta))...
                              ' norm(S,1):' num2str(norm(S,1))...
                              ' |T*S-S*T|_F2:' num2str(norm(Theta*S-S*Theta,'fro')^2)])
                        [S_hat,out] = estimate_S(C,model,regs);
                        S_hat = S_hat/max(max(S_hat));
                        Pr = out.Pr/max(max(out.Pr));
                        if isnan(S_hat)
                            S_hat = zeros(N);
                        end
                        fsc_Sg(gt,gs,ns,ct,na) = fscore(S,S_hat);
                        err_Sg(gt,gs,ns,ct,na) = norm(S - S_hat,'fro')^2/norm(S,'fro')^2;
                        err_Tg(gt,gs,ns,ct,na) = norm(Theta-Pr,'fro')^2/norm(Theta,'fro')^2;
                    end
                end 
            end 
        end
    end
    %all_spr(ng,:,:)= out.spr;
    %all_differ(ng,:,:)=out.differ;
    % figure()
    % subplot(121)
    % imagesc(out.Pr)
    % title('Theta est')
    % colorbar()
    % subplot(122)
    % imagesc(Theta)
    % title('Theta true')
    % colorbar()
    fsc_S(ng,:,:,:,:,:) = fsc_Sg;
    err_T(ng,:,:,:,:,:) = err_Tg;
    err_S(ng,:,:,:,:,:) = err_Sg;
    
end

