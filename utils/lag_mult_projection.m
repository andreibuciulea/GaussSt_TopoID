function A = lag_mult_projection(B,delta)
    nB = norm(B,"fro"); 
    if nB > delta
        A = delta/nB*B;
    else
        A = B;
    end
end