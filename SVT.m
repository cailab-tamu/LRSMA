function Y = SVT(L,lamda)
lamda = max(lamda,0);
[U,W,V] = svd(L,'econ');
VT = V';
w = diag(W);
ind = find(w>lamda);
W = diag(w(ind)-lamda);
Y = U(:,ind)*W*VT(ind,:);
end
