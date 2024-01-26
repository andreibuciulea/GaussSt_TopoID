clear 
names = {'AAPL','AMZN','GOOG','META','MSFT','NFLX','UBER'};

N = numel(names); 
for n = 1:N
    name = [names{n} '.csv'];
    data = readtable(name);
    X = table2array(data(:,6));
    X1 = X(1:end-1);
    X2 = X(2:end);
    log_returns(:,n) = log(X1./X2);
    precios(:,n) = table2array(data(:,2)); 
    open(:,n) = table2array(data(:,2)); 
    close(:,n) = table2array(data(:,6)); 
end
%%
%generate 200 graph estimations from 200 sample correlation matrix with M =
%30;
T=200;
M = 30;
N = size(log_returns,2);
C = zeros(N,N,T);
%models = {'GST-fast'};
models = {'GSR','GL','GST-fast'};
Mo = numel(models);
max_iters=10;
all_S = zeros(N,N,Mo,T);
all_Pr = zeros(N,N,Mo,T);
select_X = zeros(N,M,T);
prms = struct('sig_type',' ','g_type', ' ','C_type',' ','max_iters',max_iters,'verbose',false);
S_true = ones(N)-eye(N);
for t = 1:T
    t
    X_aux = log_returns(t:t+M-1,:);
    select_X(:,:,t) = X_aux'; 
    C_aux = cov(X_aux);
    C(:,:,t) = inv(sqrt(diag(diag(C_aux))))*C_aux*inv(sqrt(diag(diag(C_aux))));
    for m = 1:Mo
        regs = get_reg([models{m} '-finance'],prms);
        regs.S_true=S_true;
        [~,out] = estimate_S(C(:,:,t),models{m},regs,S_true);
        all_S(:,:,m,t) = out.S_hat;
        all_Pr(:,:,m,t) = out.Pr;
    end
end


%%
figure
for m = 1:3
    for t = 1:5
        subplot(3,5,(m-1)*5+t)
        imagesc(all_S(:,:,m,t))
        title(models{m})
        colorbar()
    end
end



%%
load('financial_data_4_algorithms.mat')
% figure()
% i=1;
% for k=[1,50,175,199]
%     S0 = all_S(:,:,k);
%     S0 = (S0+S0')/2;
%     S0 = S0/max(max(S0));
%     subplot(2,4,i)
%     imagesc(S0)
%     colorbar
%     subplot(2,4,i+4)
%     G = graph(S0,names);
%     LWidths = 4.5*G.Edges.Weight/max(G.Edges.Weight);
%     plot(G,'LineWidth',LWidths,'Layout','circle')
%     title(datestr(table2array(data(k+30,1))))
%     i=i+1;
% end

models = {'GSR','GL','GST-fast','SGL k comp'};
%models = {'GST-fast'};
Mo = numel(models);
alg_con_ind = zeros(Mo,T);
for t = 1:T
    for m = 1:Mo
        S = all_S(:,:,m,t);
        D = diag(sum(S));
        L = D-S;
        L = D^(-1/2)*L*D^(-1/2);
        eigL = sort(eig(L));
        alg_con_ind(m,t) = eigL(2);
    end
end
figure('Position',[100,100,1500,400])
plot(table2array(data(32:end,1)),alg_con_ind(1:3,:),LineWidth=3)
%title('Algebraic conectivity indicator \lambda_2')
legend('GSR','GL','GGSR')
grid on
ax = gca;
ax.FontSize = 25;
grid on
%
invest = 1000;
%gasto 1000 euros por empresa
total = invest/7;
%número de acciones de cada empresa
num_acc = total./precios(1,:);
%total de ganancias en funcion del precio de las acciones
ganancias1 = num_acc*precios';

figure(2)
%plot(table2array(data(:,1)),ganancias1)
%grid on

%si está por encima del umbral, entonces no compro, si está por debajo del
%umbral entonces si compro y vendo en el mismo día, entiendo.


%alg_con --> es la variable que tengo que umbralizar
TH = 1;
valors = zeros(4,numel(TH,1));
gan = zeros(Mo,T+1);
ths = [1.1,0.75,1,0.9];
%ths = ones(4,1)*0.9;

for m = 1:Mo
    alg_con = alg_con_ind(m,:);
    ii=1;
    for th = TH
        %th = 3;
        bin_algcon = alg_con > ths(m);
        ganancias2 = zeros(T+1,1);
        ganancias2(1) = invest;
        ganancias3 = ganancias2;
        for k = 1:T 
            %empiezo por el día 31
            %compro con el precio de apertura en 2 
            num_acc2 = (ganancias2(k)/7)./close(k+M-1,:);
            %vendo con el precio de cierre en 6
            ganancias2(k+1) = num_acc2*close(k+M,:)';
            if bin_algcon(k)
                ganancias3(k+1) = ganancias3(k);
            else
                %compro con el precio de apertura en 2 
                num_acc3 = (ganancias3(k)/7)./close(k+M-1,:);
                %vendo con el precio de cierre en 6
                ganancias3(k+1) = num_acc3*close(k+M,:)';
            end
        end
        valors(m,ii) = ganancias3(end);
        ii=ii+1;
    end
    gan(m,:) = ganancias3';
end


% plot(TH,valors','LineWidth',2);legend(models);
% hold on
% plot(table2array(data(31:end,1)),ganancias2)
% hold on 
% plot(table2array(data(31:end,1)),ganancias3)
% grid on
% legend('Strategy 1','Strategy 2')


figure('Position',[100,100,1500,400])

plot(table2array(data(31:end,1)),gan(1:3,:)','LineWidth',3);
hold on
plot(table2array(data(31:end,1)),ganancias2,'LineWidth',3)
legend('GSR','GL','GGSR','Strategy I')
%title(['Threshold: best'])
ax = gca;
ax.FontSize = 25;
grid on


