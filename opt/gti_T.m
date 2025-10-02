function [T,Q,Y] = gti_T(C,S,T,Q,Y,prms)
    beta = prms.beta;
    delta = prms.delta;
    L2 = prms.L2;
    In = eye(size(C,1));
    %%%%% For the second subproblem we update T, Q, and Y
    R = 1;
    for r = 1:R
        %gradient of T
        fT = S*S*T+T*S*S-2*S*T*S+Q*S-S*Q + 1/beta*(S*Y-Y*S);
        % define matrix W 
        W = T - 1/L2*(beta*fT + 1*C); %para MRF hay que poner 1*C
        %eigenvalue decomposition of W
        [U,La] = eig(W);
        %update T
        T = U*((La + sqrt(La^2+4/L2*In))/2)*U';
        %update Q 
        Qp = S*T-T*S+1/beta*Y;
        Q = lag_mult_projection(Qp,delta);
        %update Y 
        Y = Y + beta*(S*T-T*S-Q);
    end