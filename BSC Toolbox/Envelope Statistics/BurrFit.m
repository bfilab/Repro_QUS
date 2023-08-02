function [b lambda] = BurrFit(Env)
Env = Env(:);
[N, X] = hist(Env,200);
DeltaX = X(2)-X(1);

PDF_Env = N/sum(N)/DeltaX;

M_env = sum(PDF_Env.*X*DeltaX);

fun_env = @(x)CostBurrRatio(x,M_env^2/sum(PDF_Env.*(X.^2)*DeltaX));

b = fminsearch(fun_env, 3);
lambda = sqrt(sum(PDF_Env.*(X.^2)*DeltaX)*(b-2));


end
