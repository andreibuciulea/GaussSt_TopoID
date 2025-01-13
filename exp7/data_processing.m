% Load financial data and process returns
names = {'AAPL', 'AMZN', 'GOOG', 'META', 'MSFT', 'NFLX', 'UBER'};
N = numel(names);
for n = 1:N
    data = readtable([names{n} '.csv']);
    X = table2array(data(:, 6));
    log_returns(:, n) = log(X(1:end-1) ./ X(2:end));
    precios(:, n) = table2array(data(:, 2));
    close(:, n) = table2array(data(:, 6));
end

% Generate graph estimations
T = 200; % Number of graph samples
M = 30;  % Sample size for correlation estimation
Mo = 3;  % Number of models
models = {'GSR', 'GL', 'GST-fast'};
max_iters = 10;

all_S = zeros(N, N, Mo, T);
all_Pr = zeros(N, N, Mo, T);
select_X = zeros(N, M, T);

for t = 1:T
    X_aux = log_returns(t:t+M-1, :);
    select_X(:, :, t) = X_aux';
    C_aux = cov(X_aux);
    C_norm = inv(sqrt(diag(diag(C_aux)))) * C_aux * inv(sqrt(diag(diag(C_aux))));
    for m = 1:Mo
        regs = get_reg([models{m} '-finance'], struct('max_iters', max_iters, 'verbose', false));
        [~, out] = estimate_S(C_norm, models{m}, regs, ones(N) - eye(N));
        all_S(:, :, m, t) = out.S_hat;
        all_Pr(:, :, m, t) = out.Pr;
    end
end

% Plot estimated graphs
figure;
for m = 1:3
    for t = 1:5
        subplot(3, 5, (m-1)*5 + t);
        imagesc(all_S(:, :, m, t));
        title(models{m});
        colorbar();
    end
end

% Compute algebraic connectivity
alg_con_ind = zeros(Mo, T);
for t = 1:T
    for m = 1:Mo
        S = all_S(:, :, m, t);
        L = diag(sum(S)) - S; % Laplacian
        L_norm = inv(sqrt(diag(diag(L)))) * L * inv(sqrt(diag(diag(L))));
        eigL = sort(eig(L_norm));
        alg_con_ind(m, t) = eigL(2);
    end
end

% Plot algebraic connectivity
figure('Position', [100, 100, 1500, 400]);
plot(1:T, alg_con_ind', 'LineWidth', 3);
legend(models);
grid on;

% Simulate investment strategies
invest = 1000; % Initial investment
TH = [1.1, 0.75, 1, 0.9];
gan = zeros(Mo, T+1);
for m = 1:Mo
    alg_con = alg_con_ind(m, :);
    bin_algcon = alg_con > TH(m);
    ganancias = zeros(T+1, 1);
    ganancias(1) = invest;
    for k = 1:T
        if ~bin_algcon(k)
            num_acc = (ganancias(k) / 7) ./ close(k+M-1, :);
            ganancias(k+1) = num_acc * close(k+M, :)';
        else
            ganancias(k+1) = ganancias(k);
        end
    end
    gan(m, :) = ganancias;
end

% Plot strategy results
figure('Position', [100, 100, 1500, 400]);
plot(1:T+1, gan', 'LineWidth', 3);
legend('GSR', 'GL', 'GST-fast');
grid on;
