function A = generate_RBF_graph(N,T,s, connected)
max_tries = 100;
    A = gaussian_graph(N,T,s);
    if connected %generate connected graph
        k = 1;
        lambdas = eig(diag(sum(A))-A);
        while((abs(lambdas(2)) < 1e-6) && (k < max_tries))
            A = gaussian_graph(N,T,s);
            k = k+1;
            lambdas = eig(diag(sum(A))-A);
        end
        if k >= max_tries
            disp('ERR: Generated unconnected RBF graph')
        end
    end
end