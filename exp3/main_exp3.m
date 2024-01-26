clear all
clc
addpath('../utils');
addpath('../opt');
addpath(genpath('../../global_utils'));

g_type = {'BA','SW','SBM'};
%g_type = {'SW'};
g_size = [20];
n_samples = [1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6];
%n_samples = [1e2, 1e3, 1e4, 1e5, 1e6];
%n_samples = [1e2,1e3,1e4,1e5,1e6];
%n_samples = [1e6];
C_types = {'SSEM','Poly'};
%C_types = {'Poly'};
%models = {'GGSR','GGSR-jiaxi-v1'};
models = {'GGSR-jiaxi-v1'};

rng(2)
nG = 128;
Gt = numel(g_type);
Gs = numel(g_size);
nS = numel(n_samples);
Ct = numel(C_types);
nA = numel(models);

%%%Graph gen params
N = 20;
ER_p = 0.1;
RBF_T = 0.9;RBF_s = 0.5;RBF_conn = true;
BA_m = 2;
SW_K = 2;SW_Beta=0.15;
SBM_k = 4;SBM_p = 0.8;SBM_q = 0.05;
norm_L = true;L_bin = true;

%%%Signal gen params 
max_iters = 20;
verbose = false;
norm_noise = true;
sampled = true;
L = 3;
sigma = 0;
    
fsc = zeros(nG,Gt,Gs,nS,Ct,nA);
fronorm = zeros(nG,Gt,Gs,nS,Ct,nA);
tic
parfor ng = 1:nG %number of graphs
    disp(['Graph: ' num2str(ng)])
    %initialize variables
    fsc_g = zeros(Gt,Gs,nS,Ct,nA);
    fronorm_g = zeros(Gt,Gs,nS,Ct,nA);
    %general parameters
    prms = struct('verbose',verbose,'norm_noise',norm_noise,'sampled',sampled,...
              'max_iters',max_iters,'sigma',sigma,'L',L);
    g_prms = struct('N',N,'ER_p',ER_p,'RBF_T',RBF_T,'RBF_s',RBF_s,'RBF_conn',RBF_conn,...
    'BA_m',BA_m,'SW_K',SW_K,'SW_Beta',SW_Beta,'SBM_k',SBM_k,'SBM_p',SBM_p,'SBM_q',SBM_q,...
    'norm_L',norm_L,'L_bin',L_bin);    
    for gt = 1:Gt %number of graph types
        g_prms.g_type = g_type{gt}; 
        prms.g_type = g_type{gt};
        for gs = 1:Gs %number of graph sizes
            g_prms.N = g_size(gs);
            S = generate_graph(g_prms).A;
            for ns = 1:nS %number of samples
                prms.M = n_samples(ns);
                for ct = 1:Ct %number of types of C
                    prms.sig_type = C_types{ct};
                    gsout = generate_graph_signals(S, prms);
                    C = gsout.C;
                    %C = C/max(max(C));
                    for na = 1:nA %number of algorithms
                        model = models{na};
                        regs = get_reg(model,prms);
                        regs.S_true=S;
                        [S_hat,~] = estimate_S(C,model,regs);
                        S_hat = S_hat/max(max(S_hat));
                        fsc_g(gt,gs,ns,ct,na) = fscore(S,S_hat);
                        fronorm_g(gt,gs,ns,ct,na) = norm(S-S_hat,'fro')^2/norm(S,'fro')^2;
                    end
                end 
            end 
        end
    end
    fsc(ng,:,:,:,:,:) = fsc_g;
    fronorm(ng,:,:,:,:,:) = fronorm_g;
end
toc
%% Figure 4 
load('data_exp3_v2.mat')
%%
mrkt = {'+-','+--','x-','x--','*-','*--'};

tp = 'mean';
if strcmp(tp,'mean')
    %aux1 = squeeze(nanmean(fsc));
    aux1 = squeeze(mean(fronorm));
else
    aux1 = squeeze(sum(fsc==1)/nG);
end

figure('Position',[100,100,800,550])
i = 1;
for k = [2,3,1] %1:numel(g_type)
    loglog(n_samples,aux1(k,:,2),mrkt{k*2-1},'LineWidth',4,'MarkerSize',14)
    lgd{i*2-1} = ['GGSR ' g_type{k} ' Poly'];
    hold on
    loglog(n_samples,aux1(k,:,1),mrkt{k*2},'LineWidth',4,'MarkerSize',14)
    lgd{i*2} = ['GGSR ' g_type{k} ' SSEM'];
    hold on
    i = i+1;
end
%si es grafo no es muy esparse esto funciona peor
lg = legend(lgd,'FontSize',18,'FontWeight','bold', 'Interpreter','latex')
set(lg,'color','none');
xlabel('Number of samples','FontSize',28,'FontWeight','bold', 'Interpreter','latex')
if strcmp(tp,'mean')
    ylabel('Normalized mean error','FontSize',28,'FontWeight','bold', 'Interpreter','latex')
else
    ylabel('Ratio of recovered graphs','FontSize',28,'FontWeight','bold', 'Interpreter','latex')
end
ax = gca;
ax.FontSize = 24;
%title('est-S-GSTRw','FontWeight','bold','FontSize',12)
grid on

xticklabels({'10^2','10^3','10^4','10^5','10^6'})
xticks([1e2,1e3,1e4,1e5,1e6])

yticklabels({'10^{-4}','10^{-3}','10^{-2}','10^{-1}','10^{0}'})
yticks([1e-4,1e-3,1e-2,1e-1,1e0,])