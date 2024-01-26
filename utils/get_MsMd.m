function [Ms,Md] = get_MsMd(N)
    R = N*(N+1)/2; 
    W = tril(ones(N))';
    idx = find(W(:) ~= 0);
    Ms = zeros(R,N^2);
    for k = 1:R
        Ms(k,idx(k)) = 1;
    end
    Md = zeros(N^2,R);
    k = 1;
    for i = 1:N
        for j = 1:i
            %disp([num2str(j) num2str(i)]);
            id1 = N*(i-1)+j;
            id2 = N*(j-1)+i;
            Md(id1,k) = 1;
            Md(id2,k) = 1;
            k = k+1;
        end
    end

end