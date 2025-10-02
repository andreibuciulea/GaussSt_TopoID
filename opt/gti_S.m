function [S,Z,P] = gti_S(T,S,Z,P,prms)  
    beta = prms.beta;
    eta = prms.eta;
    rho = prms.rho;
    delta = prms.delta;
    L1 = prms.L1;
    constraint = prms.constraint;
    N = size(T,1);
    In = eye(N);
    R = prms.mit; 
    %%%%% For the first subproblem we update S, P, and Z
    la_max = 1*max(abs(eig(T)));
    L1 = (4*beta*la_max^2);
    if sum(sum(S)) == 0
        S = generate_connected_ER(N,0.1);
    end
    for r = 1:R
        %gradient of S
        gS = T*T*S + S*T*T - 2*T*S*T + P*T - T*P + 1/beta*(T*Z - Z*T);
        %update S before projection
        Sp = S - 1/L1*(rho*In + eta*S + beta*gS); 
        %projection for S
        %S = projection_symmetry_simplex(Sp,constraint);
        Sp(Sp<0) = 0;
        S = (Sp+Sp')/2;
        S = S/max(max(S));
        S = S-diag(diag(S));
        %figure(2);imagesc(S);colorbar();
        %S = S/max(max(S));
        %update for P
        Pp = T*S-S*T+1/beta*Z;
        P = lag_mult_projection(Pp,delta);
        %update for Z
        Z = Z + beta*(T*S-S*T-P);
    end

end