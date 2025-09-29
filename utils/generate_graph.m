function out = generate_graph(pr) 
    g_type = pr.g_type;
    N = pr.N;
    if isfield(pr,'norm_L')
        norm_L = pr.norm_L;
    else
        norm_L = true;
    end
    if isfield(pr,'L_bin')
        L_bin = pr.L_bin;
    else
        L_bin = true;
    end
    
    A = [];
    L = [];
    switch g_type
        case 'ER'
            A = generate_connected_ER(N,pr.ER_p);
        case 'TREE'
            A = preferential_attachment_graph(N,1);
        case 'BA'
            A = preferential_attachment_graph(N,pr.BA_m);
        case 'RBF'
            A = generate_RBF_graph(N,pr.RBF_T,pr.RBF_s,pr.RBF_conn);
        case 'RING'
            A = small_world(N,pr.RING_K,0);
        case 'SW'
            %K = pr.K;beta = pr.beta;A = small_world(N,K,beta);
            h2 = WattsStrogatz(N,pr.SW_K,pr.SW_Beta);
            A = full(adjacency(h2));
        case 'SBM'
            prms = struct('p', pr.SBM_p,'q',pr.SBM_q);
            G = gsp_stochastic_block_graph(N, pr.SBM_k,prms);
            A = full(G.A);
        otherwise
            disp('ERR: Unkown graph type')
            return
    end
    while min(abs(eig(eye(N)-A))) <= 1e-1  
        A = generate_graph(pr).A;%for SSEM signals
    end
    if L_bin
        A = double(A > 0.5);
    end
    L = diag(sum(A,2)) - A;
    if norm_L
        L = L/trace(L)*N;
        %L = L/norm(L,'fro');
    else
        % sigma = sigma*sqrt(N/trace(L));
    end
    out.A = A;
    out.L = L;
end