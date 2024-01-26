function F = objective_function(Theta, S, Sigma, rho, beta, Y)
[A, FLAG]=chol(Theta,'lower'); % FLAG = 1 if Theta is not PSD
lg_det =2*sum(log(diag(A)));
fun = -1*lg_det + dot(Theta(:), Sigma(:)) + rho*sum(sum(abs(S))) + 0.5*beta*norm(Theta*S-S*Theta,'fro')^2 + trace(Y*(Theta*S-S*Theta));
F.value = fun;
F.flag = FLAG;
end

