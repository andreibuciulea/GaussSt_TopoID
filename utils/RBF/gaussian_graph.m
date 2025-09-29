function G = gaussian_graph(N,T,s)
    XCoords = rand(N,1);
    YCoords = rand(N,1);
    d = distanz(XCoords,YCoords);
    W = exp(-d.^2/(2*s^2)); 
    W(W<T) = 0; % Thresholding to have sparse matrix
    W = 0.5*(W+W');
    G = W-diag(diag(W));
end
