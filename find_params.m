%
addpath('utils');
addpath('global_utils')
addpath('opt');
rng(2);

N = 20;
M = 1e6;
g_type = 'ER';

%%%%%%%%%%%% Signal generation parameters
sig_type = 'MRF';
L = 3;
sampled = true;
max_iters = 50;
prms = struct('sig_type',sig_type,'L',L,'M',M,'sampled',sampled,...
              'max_iters',max_iters);
%%%%%%%%%%%%
%%%Graph gen params
    norm_L = true;
    L_bin = true;
    %ER
    p = 0.1/(N/20); %grafos igual de esparsos
    %BA
    m=2; %con m=1 es demasiado esparso y no funciona bien
    %RBF
    connected=true;
g_prms = struct('norm_L',norm_L,'L_bin',L_bin,'p',p,'m',m,'connected',connected,'N',N,'g_type',g_type); 


nG = 64;
model = 'GGSR-j';
max_iters = 50;


A0 = 1; %rho  %probar valores de rho
A1 = 1e-4; %delta
A2 = 1e-5; %beta 1e-5;
A3 = 1e-4; %eta
A4 = 1e-2; %L1
A5 = 0.2; %L2 %probar otros valores de L2 entre 0.1 y 2

n_a = numel(A0);n_b = numel(A1);n_c = numel(A2);n_d = numel(A3);n_e = numel(A4);n_f = numel(A5);
n_tot = n_a*n_b*n_c*n_d*n_e*n_f;
all_a = cell(1,n_tot);all_b=cell(1,n_tot);all_c=cell(1,n_tot);
all_d = cell(1,n_tot);all_e=cell(1,n_tot);all_f=cell(1,n_tot);
i = 1;
for a = 1:n_a
    for b = 1:n_b
        for c = 1:n_c
            for d = 1:n_d
                for e = 1:n_e
                    for f = 1:n_f
                        all_a{i} = A0(a);
                        all_b{i} = A1(b);
                        all_c{i} = A2(c);
                        all_d{i} = A3(d);
                        all_e{i} = A4(e);
                        all_f{i} = A5(f);
                        i = i+1;
                    end
                end
            end
        end
    end
end

regs = struct('rho',all_a,'delta',all_b,'beta',all_c,'eta',all_d,'L1',all_e,'L2',all_f,'max_iters',max_iters);
f_score = zeros(nG,2,n_tot);
errS = zeros(nG,n_tot);
errT = zeros(nG,n_tot);

tic

parfor g=1:nG
    disp(['Graph: ' num2str(g)])
    S = generate_graph(g_prms).A;
    out = generate_graph_signals(S, prms);
    C = out.C;
    Theta = out.C_inv;
    Theta = Theta/max(max(Theta));
    C = C/max(abs(eig(C)));

    fsc_g = zeros(n_tot,1);fsc_gT = zeros(n_tot,1);
    err_g = zeros(n_tot,1);err_gT = zeros(n_tot,1);
    
    for t = 1:n_tot
        [S_hat,out] = estimate_S(C,model,regs(t));
        Pr = out.Pr; Pr = Pr/max(max(Pr));
        %fsc_g(t) = fscore(S,S_hat);
        err_g(t) = norm(S_hat - S,'fro')^2/norm(S,'fro')^2;
        err_gT(t) = norm(Pr - Theta,'fro')^2/norm(Theta,'fro')^2;      
        %fsc_gT(t) = fscore(Theta,Pr);
    end
    %f_score(g,1,:) = fsc_g;
    %f_score(g,2,:) = fsc_gT;
    errS(g,:) = err_g;
    errT(g,:) = err_gT;

end
t = toc;
disp(['--- ' num2str(t/3600) ' hours'])

%% parfor on regularizers
for g=1:nG
    disp(['Graph: ' num2str(g)])
    S = generate_graph(g_prms).A;
    out = generate_graph_signals(S, prms);
    C = out.C;
    Theta = out.C_inv;
    Theta = Theta/max(max(Theta));
    C = C/max(abs(eig(C)));

    fsc_g = zeros(n_tot,1);fsc_gT = zeros(n_tot,1);
    err_g = zeros(n_tot,1);err_gT = zeros(n_tot,1);
    
    parfor t = 1:n_tot
        [S_hat,out] = estimate_S(C,model,regs(t));
        Pr = out.Pr; Pr = Pr/max(max(Pr));
        fsc_g(t) = fscore(S,S_hat);
        err_g(t) = norm(S_hat - S,'fro')^2/norm(S,'fro')^2;
        err_gT(t) = norm(Pr - Theta,'fro')^2/norm(Theta,'fro')^2;      
        fsc_gT(t) = fscore(Theta,Pr);
    end
    f_score(g,1,:) = fsc_g;
    f_score(g,2,:) = fsc_gT;
    err(g,1,:) = err_g;
    err(g,2,:) = err_gT;

end
t = toc;
disp(['--- ' num2str(t/3600) ' hours'])

%%
figure()
aux1 = squeeze(sum(fsc==1));
%aux1 = squeeze(mean(fsc));
imagesc(aux1)
colorbar()
set(gca,'XTick',[1:numel(lambdas1)],'XTickLabel', lambdas1);
set(gca,'YTick',[1:numel(rhos)],'YTickLabel', rhos);
set(gca,'XTickLabelRotation',45)
xlabel('lambdas1')
ylabel('rhos')
title('est S GST M = 1e6')

%%
[X,Y] = meshgrid(log10(lambdas),log10(rhos));
surf(X,Y,squeeze(sum(fsc==1)))
%%
clear 
close all
addpath('utils');
addpath('opt');
addpath('Graphical Lasso');
rng(1);


%%%%%%%%%%%% Signal generation parameters
N = 20;
M = 1e6;
g_type = 'ER';
C_type = 'Poly';
L = 5;
sigma = 0;
norm_noise = true;
sampled = true;
max_iters = 10;
verbose = false;
prms = struct('C_type',C_type,'g_type',g_type,'L',L,'verbose',verbose,...
              'sigma',sigma,'M',M,'norm_noise',norm_noise,'sampled',sampled,...
              'max_iters',max_iters);
%%%%%%%%%%%%
%%%Graph gen params
norm_L = true;
L_bin = true;
%ER
p=0.1;
%BA
m=2; %con m=1 es demasiado esparso y no funciona bien
%RBF
connected=true;
g_prms = struct('norm_L',norm_L,'L_bin',L_bin,'p',p,'m',m,'connected',connected,'N',N,'g_type',g_type); 

%%%%%%%%%%%%
%%%Generate graphs
nG = 32;
all_S = cell(1,nG);all_C =cell(1,nG);
for g = 1:nG
    S = generate_graph(g_prms); 
    [~,~,all_C{g},~] = generate_graph_signals(S, prms);
    all_S{g} = S;
end
model = 'GST-fast';

kappa= 1e-2;
tol=5e-2;
alpha=1e-4;
rhos = 8;
lambdas0 = 0;
lambdas1 = 1e-6;%5
lambdas2 = 1e-3;%5
lambdas3 = 0;
verbose = 0;
is_MRF = 0;
%generate all regs combinations
n_ro = numel(rhos);n_la0 = numel(lambdas0);n_la1 = numel(lambdas1);n_la2 = numel(lambdas2);n_la3 = numel(lambdas3);
n_tot = n_ro*n_la0*n_la1*n_la2*n_la3;
all_ro = cell(1,n_tot);all_la0 = cell(1,n_tot);all_la1 = cell(1,n_tot);all_la2 = cell(1,n_tot);all_la3 = cell(1,n_tot);
all_verb = num2cell(logical(ones(1,n_tot)*verbose));
all_iters = num2cell(ones(1,n_tot)*max_iters);
all_isMRF = num2cell(logical(ones(1,n_tot)*is_MRF));

i = 1;
for r = 1:n_ro
    for la0 = 1:n_la0
        for la1 = 1:n_la1
            for la2 = 1:n_la2
                for la3 = 1:n_la3
                    all_ro{i} = rhos(r);
                    all_la0{i} = lambdas0(la0);
                    all_la1{i} = lambdas1(la1);
                    all_la2{i} = lambdas2(la2);
                    all_la3{i} = lambdas3(la3);
                    i = i+1;
                end
            end
        end
    end
end

errS = zeros(nG,n_tot);
fscS = zeros(nG,n_tot);
errT = zeros(nG,n_tot);
fscT = zeros(nG,n_tot);
regs = struct('rho',all_ro,'la0',all_la0,'la1',all_la1,'la2',all_la2,'la3',all_la3,...
                'verbose',all_verb,'kappa',kappa,'tol',tol,'max_iters',all_iters,...
                'alpha',alpha,'is_MRF',all_isMRF);
tic
parfor ng = 1:nG
    disp(num2str(ng))
    S = generate_graph(g_prms);
    [~,~,C,Theta] = generate_graph_signals(S, prms);
    C = C/max(abs(eig(C)));
    %S = all_S{ng};C = all_C{ng};
    err_tS = zeros(1,n_tot);
    fsc_tS = zeros(1,n_tot);
    err_tT = zeros(1,n_tot);
    fsc_tT = zeros(1,n_tot);
    for t=1:n_tot
        [S_hat,out] = estimate_S(C,model,regs(t),S);  
        fsc_tS(t) = fscore(S,S_hat);
        err_tS(t) = norm(out.S_hat/norm(out.S_hat,'fro') - S/norm(S,'fro'),"fro");
        err_tT(t) = norm(out.Pr/norm(out.Pr,'fro') - Theta/norm(Theta,'fro'),'fro');
        fsc_tT(t) = fscore(Theta,out.Pr);
    end
    errS(ng,:) = err_tS
    fscS(ng,:) = fsc_tS;
    errT(ng,:) = err_tT
    fscT(ng,:) = fsc_tT;
    figure()
    subplot(121)
    imagesc(out.Pr)
    title('Theta est')
    colorbar()
    subplot(122)
    imagesc(Theta)
    title('Theta true')
    colorbar()
end
t = toc;
disp(['--- ' num2str(t/3600) ' hours'])


%%

figure()
plot(errS)
hold on
plot(errT)

