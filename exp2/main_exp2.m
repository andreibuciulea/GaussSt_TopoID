addpath(genpath('../utils'));
addpath(genpath('../opt'));
addpath(genpath('../../global_utils'));

rng(6);

g_type = {'ER'};
g_size = [20];  
%n_samples = [1e2, 3e2, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]; 
n_samples = [1e2, 1e3, 1e4, 1e5, 1e6]; 
C_types = {'Poly'};
%sigma_values = [0.05,0.1,0.2,0.3];
sigma_values = [0.05, 0.2];
models = {'GL','GSR','GGSR','GGSR-jiaxi-v1'};

n_samples = [1e6];
models = {'GSR','GGSR','GGSR-jiaxi-v1'};
nG = 64;
Gt = numel(g_type);
Gs = numel(g_size);
nS = numel(n_samples);
Ct = numel(C_types);
Sv = numel(sigma_values);
nA = numel(models);

%%%Graph gen params
ER_p = 0.1;
RBF_T = 0.9;RBF_s = 0.5;RBF_conn = true;
BA_m = 2;
SW_K = 2;SW_Beta=0.15;
SBM_k = 4;SBM_p = 0.8;SBM_q = 0.05;
norm_L = true;L_bin = true;

%%%Signal gen params 
    max_iters = 100;
    verbose = false;
    norm_noise = true;
    sampled = true;
    sigma = 0;
    L = 5;
    
fsc = zeros(nG,Gt,Gs,nS,Ct,Sv,nA);
fronorm = zeros(nG,Gt,Gs,nS,Ct,Sv,nA);
fr = zeros(nG,Gt,Gs,nS,Ct,Sv);
tic
parfor ng = 1:nG %number of graphs
    disp(['Graph: ' num2str(ng)])
    %initialize variables
    fsc_g = zeros(Gt,Gs,nS,Ct,Sv,nA);
    fronorm_g = zeros(Gt,Gs,nS,Ct,Sv,nA);
    fr_g = zeros(Gt,Gs,nS,Ct,Sv);
    %general parameters
    prms = struct('verbose',verbose,'norm_noise',norm_noise,'sampled',sampled,...
              'max_iters',max_iters,'sigma',sigma,'L',L);
    g_prms = struct('ER_p',ER_p,'RBF_T',RBF_T,'RBF_s',RBF_s,'RBF_conn',RBF_conn,...
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
                    for sv = 1:Sv 
                        prms.sig_type = C_types{ct};
                        prms.sigma = sigma_values(sv);
                        gsout = generate_graph_signals(S, prms);
                        C = gsout.C;
                        for na = 1:nA %number of algorithms
                            model = models{na};
                            regs = get_reg(model,prms);
                            regs.S_true = S;
                            [S_hat,out] = estimate_S(C,model,regs);
                            fsc_g(gt,gs,ns,ct,sv,na) = fscore(S,S_hat);
                            fronorm_g(gt,gs,ns,ct,sv,na) = norm(S-S_hat,'fro')^2/norm(S,'fro')^2;

                            %%% Precision matrix for GST
                            if strcmp(model,'est-S-GST')
                                SN = norm(S,'fro');
                                S_pr = out.Pr-diag(diag(out.Pr));
                                fr_g(gt,gs,ns,ct,sv) = norm(S-S_pr,'fro')^2/(SN)^2;
                            end
                        end
                    end
                end 
            end 
        end
    end
    fsc(ng,:,:,:,:,:,:) = fsc_g;
    fronorm(ng,:,:,:,:,:,:) = fronorm_g;
    fr(ng,:,:,:,:,:) = fr_g;
end
toc
%% Figure 5
load('data_exp2_v3.mat')
%%
%aux1 = squeeze(mean(fsc));
%aux1 = squeeze(sum(fsc==1)/nG);
aux1 = squeeze(mean(fronorm));
mrkt = {'^-','^--','s-','s--','*-','*--'};

c = 1;
a = 1;
b = 3;
models= {'GL','GSR','GGSR'};
figure('Position',[100,100,800,550])
for k = 1:numel(models)
    loglog(n_samples,squeeze(aux1(:,c,a,k)),mrkt{2*k-1},'LineWidth',3,'MarkerSize',14)
    lgd{2*k-1} = [models{k} ' ($\sigma$=' num2str(sigma_values(a)) ')'];
    hold on
    loglog(n_samples,squeeze(aux1(:,c,b,k)),mrkt{2*k},'LineWidth',3,'MarkerSize',14)
    lgd{2*k} = [models{k} ' ($\sigma$=' num2str(sigma_values(b)) ')'];
    hold on
end
%%%%%%%%%%%%%%%
% aux2 = squeeze(mean(fr));
% loglog(n_samples,squeeze(aux2(:,c,a)),'o-.','LineWidth',3,'MarkerSize',14)
% lgd{7} = ['GGSR-{\boldmath$\Theta$} ($\sigma$=' num2str(sigma_values(a)) ')'];
% aux2 = squeeze(mean(fr));
% loglog(n_samples,squeeze(aux2(:,c,b)),'o-.','LineWidth',3,'MarkerSize',14)
% lgd{8} = ['GGSR-{\boldmath$\Theta$} ($\sigma$=' num2str(sigma_values(b)) ')'];

%%%%%%%%%%%%%%%
lg = legend(lgd,'FontSize',15,'FontWeight','bold', 'Interpreter','latex')
set(lg,'color','none');
xlabel('Number of samples','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
ylabel('Normalized mean error','FontSize',24,'FontWeight','bold', 'Interpreter','latex')
ax = gca;
ax.FontSize = 20;
grid on
%title(C_types{c})

c = 2;
aux1(:,c,1,1) = 1;
figure('Position',[100,100,800,550])
for k = 1:numel(models)
    loglog(n_samples,squeeze(aux1(:,c,b,k)),mrkt{2*k-1},'LineWidth',4,'MarkerSize',14)
    lgd{2*k-1} = [models{k} ' ($\sigma$=' num2str(sigma_values(b)) ')'];
    hold on
    loglog(n_samples,squeeze(aux1(:,c,a,k)),mrkt{2*k},'LineWidth',4,'MarkerSize',14)
    lgd{2*k} = [models{k} ' ($\sigma$=' num2str(sigma_values(a)) ')'];
    hold on

end
%%%%%%%%%%%%%%%%%%%%%%%%%%
% loglog(n_samples,squeeze(aux2(:,c,a)),'o-.','LineWidth',3,'MarkerSize',14)
% lgd{7} = ['GGSR-{\boldmath$\Theta$} ($\sigma$=' num2str(sigma_values(a)) ')'];
% aux2 = squeeze(mean(fr));
% loglog(n_samples,squeeze(aux2(:,c,b)),'o-.','LineWidth',3,'MarkerSize',14)
% lgd{8} = ['GGSR-{\boldmath$\Theta$} ($\sigma$=' num2str(sigma_values(b)) ')'];
%%%%%%%%%%%%%%%%%%%%%%%%
lg = legend(lgd,'FontSize',18,'FontWeight','bold', 'Interpreter','latex')
set(lg,'color','none');
xlabel('Number of samples','FontSize',28,'FontWeight','bold', 'Interpreter','latex')
ylabel('Normalized mean error','FontSize',28,'FontWeight','bold', 'Interpreter','latex')
ax = gca;
ax.FontSize = 24;
grid on
ylim([0.1 1.1])
%title(C_types{c})

yticks([0.1 0.2 0.4 0.6 0.8 1])
yticklabels({'0.1','0.2','0.4','0.6','0.8','1'})

xticklabels({'10^2','10^3','10^4','10^5','10^6'})
xticks([1e2,1e3,1e4,1e5,1e6])

