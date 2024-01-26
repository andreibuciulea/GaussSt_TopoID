function val = obj_fun(Sigma,Theta,S,rho,beta)    
   
    val = -log(det(Theta)) + trace(Sigma*Theta) + rho*norm(S,1)...
        + beta/2*norm(Theta*S-S*Theta,'fro')^2;

end