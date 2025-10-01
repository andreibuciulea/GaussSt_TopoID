function out = generate_graph_signals(S,prms)

    N = size(S,1);
    if isfield(prms,'M'); M = prms.M; 
    else; M = 1e2; end
    if isfield(prms,'L'); L = prms.L; 
    else; L = 3; end
    if isfield(prms,'sigma'); sigma = prms.sigma; 
    else; sigma = 0; end
    if isfield(prms,'sig_type'); sig_type = prms.sig_type; 
    else; sig_type = 'smooth';end
    if isfield(prms,'verbose');verbose = prms.verbose; 
    else; verbose = false; end
    if isfield(prms,'norm_noise'); norm_noise = prms.norm_noise; 
    else; norm_noise = true;end
    if isfield(prms,'sampled');sampled = prms.sampled; 
    else; sampled = false; end

    X = zeros(N,M); X_0 = X;
    C = zeros(N,N); snr = 0;

    switch sig_type  %create separate functions for each case
        case 'Poly'    
            h1 = rand(L,1); % Draw the coefficients of the first polynomial
            H1 = zeros(N,N);
            for ii = 1:L
                H1 = H1 + h1(ii)*S^(ii-1);
            end
            C = H1^2;
            C_inv = inv(C);
            %C_inv = inv(C/max(abs(eig(C))));
            %C_inv = N*C_inv/trace(C*C_inv);
            %delete next 2 lines
            %C_inv = H1^2;
            %C = inv(C_inv);
            X_iid=randn(N,M);
            X_0 = sqrtm(C)*X_iid;
        case 'MRF'
            [~,D] = eig(S);
            C_inv = (0.01 - min(diag(D)))*eye(N,N) + (0.9 + 0.1*rand(1,1))*S;
            C = inv(C_inv);
            [~,Dc] = eig(C);
            if min(diag(Dc)) < 0
                disp('no es def pos')
                [~,~,C,~] = generate_graph_signals(S,prms);
            end
            X_iid=randn(N,M);
            X_0 = sqrtm(C)*X_iid;
        case 'smooth'
            Lap = diag(sum(S,2)) - S;
            h_mu = zeros(N,1);
            alpha = 10;%alpha indicates how smooth the signal is
            [V, Lambda] = eig(Lap); 
            Lambda_inv = inv(alpha*Lambda+eye(N));
            H = mvnrnd(h_mu, Lambda_inv, M)';
            X_0 = V*H;
            C = pinv(Lap);
            C_inv = Lap;
        case 'SSEM'
            P = randn(N,M);
            X_0 = inv(eye(N)-S)*P;
            C_inv = (eye(N)-S)^2;
            C = inv(C_inv);
        otherwise
            error('ERR: Unkown signal type')
    end
    
    p_x = norm(X_0, 'fro')^2/M;
    
    % Set sigma for normalizing the noise power
    if  norm_noise
        sigma = sqrt(sigma*p_x/N);
    end
    noise = randn(N, M)*sigma;
    
    X = X_0 + noise;
    if sampled
       C = X*X'/M;
    end
    p_n = norm(noise, 'fro')^2/M;
    snr = p_x/p_n;
    if verbose
        eq_sigma = sigma^2*N/p_x;
        disp(['Mean Np: ' num2str(p_n) '   norm sigma: ' num2str(eq_sigma)])
        disp(['SNR(nat|dB): ', num2str(snr) ' | ' num2str(10*log10(snr))])
        %disp(['Mean smoothness: ' num2str(trace(X_0'*L*X_0)/M)])
    end

    out.X = X;
    out.X_0 = X_0;
    out.C = C;
    out.C_inv = C_inv;
    out.snr = snr;
end
