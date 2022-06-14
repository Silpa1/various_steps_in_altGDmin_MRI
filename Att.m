function x = Att(y)
global S2 n1 n2 q  mk m n 

w = zeros(n,q);
y_mat=reshape(y, [m,q]);
for k=1:1:q
    w(S2(1:mk(k),k),k)=y_mat(1:mk(k),k);
    x(:,k) = n*reshape( ifft2( reshape(w(:,k), [n1 n2]) ), [n,1]) ;
end
%Y1vec=reshape(X1vec,[n*q,1]);