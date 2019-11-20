clc,clear;


n = 2048;
A = -100 + 200*rand(n,n);
x_known = ones(n,1);
b = A*x_known;
P = 1:n;
u_i = 1:n;
u_i = u_i(:);
LU = A;
tic
for i = 1:n-1
     [m,p] = max(abs(LU(u_i(i:end),i)));
     u_i([p+i-1,i]) = u_i([i,p+i-1]);
     k = u_i(i+1:n);
     l = LU(k,i)/LU(u_i(i),i);
     LU(k,i:n) =  LU(k,i:n)- l.*LU(u_i(i),i:n);
     LU(k,i) = l;
end
LU = LU(u_i,:);
U = triu(LU);
L = tril(LU,-1);
L(1:n+1:end) = ones(n,1);

b_LU = b(u_i);
x = zeros(n,1);
y = zeros(n,1);

for i=1:n
    y(i) = b_LU(i) - dot(L(i,:),y);
end

x(n) = y(n)/U(n,n);
for i=n-1:-1:1
    x(i) = (y(i) - dot(U(i,i+1:end),x(i+1:end)))/U(i,i);
end
 toc
normi1 = norm(x-x_known, inf)
normi2 = norm(A*x-b,inf)
