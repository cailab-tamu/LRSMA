function [L,S] = LRSMA(D,beta1,beta2,beta3,tol)

if nargin < 5
    tol = 1e-4;
end

L = zeros(size(D));
S = L;
Z = L;
Y = L;
G = zeros(size(D,1)-1,size(D,2)); 
rho = norm(D(:))/10;
maxIter = 500;

for iter = 1:maxIter
    % Find low rank component using singlar value thresholding
    L = SVT((D-S+rho*(Z-Y/rho))/(1+rho),beta1/(1+rho));
    % Find smooth signals using flsa
    Zold = Z;
    for j = 1:size(Z,2)
        k = L(:,j)+Y(:,j)/(rho+eps);
        g = G(:,j); % subgradient of fused term
        [Z(:,j),G(:,j)] = flsa(k,g,0,beta2/(rho+eps),length(k),100,1e-10,1,6);
    end
    % Find sparse signals
    S = D - L;
    S = sign(S) .* max( abs(S)-beta3, 0);
    % Find Y
    Y = Y + rho*(L-Z);
    % Check convergence
    Res1 = norm(L(:)-Z(:));
    Res2 = rho*norm(Z(:)-Zold(:));
    th = norm(D(:))*tol;
    fprintf('Iter %d: L-Z = %f, Z-Residual = %f, Threshold = %f \n',iter,Res1,Res2,th);
    if  Res1 < th && Res2 < th
        break
    end
    
    if Res1>10*Res2
        rho = 2*rho;
    elseif Res2>10*Res1
        rho = rho/2;
    else
    end
end

end

