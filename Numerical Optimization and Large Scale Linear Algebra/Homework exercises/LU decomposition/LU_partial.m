clc,clear;

tic
n = 2048;
x_known = ones(n,1);
normi1 = zeros(10,1);
normi2 = zeros(10,1);
for i = 1:10
    %A = -100 + 200*rand(n,n);
     A = eye(n);
     idx = tril(true(size(A)), -1);
     A(idx) = -1;
     A(:,n)=1;
    b = A*x_known;
    [p, L, U] = partialLU(A);
    x = solveLU(p,L,U,b);

    normi1(i) = norm(x-x_known, inf);
    normi2(i) = norm(A*x-b,inf);
end
meanError = mean(normi1);
meanCorrection = mean(normi2);
toc
disp(['The average norm of x error is ',num2str(meanError),' and the average norm of Ax-b correction is ',num2str(meanCorrection)]);


function [p, L, U] = partialLU(A)
% partialLU  returns L, U
% matrices and the p vector after partial
% LU decomposition of A
LU = A;
n = length(A);
p = 1:n;
p = p(:);
for i = 1:n-1
    [~,pos] = max(abs(LU(i:end,i)));
    p([i, pos+i-1]) = p([pos+i-1, i]); % interchange 2 elements
    LU([pos+i-1, i], :) = LU([i, pos+i-1], :); % interchange 2 rows
    k = i+1:n;
    l = LU(k,i)/LU(i,i);
    LU(k,i:n) =  LU(k,i:n)- l.*LU(i,i:n); % compute upper part
    LU(k,i) = l; % save lower part in the same matrix
end
% seperate the L and U matrices
U = triu(LU);
L = tril(LU,-1);
L(1:n+1:end) = ones(n,1);
end

function x = solveLU(p,L,U,b)
% SolveLU returns the
% solution x of the system
% LUx = b
n=length(b);
b_LU = b(p);
x = zeros(n,1);
y = zeros(n,1);
% front substitution
for i=1:n
   y(i) = b_LU(i) - dot(L(i,:),y);
end
% back substitution
x(n) = y(n)/U(n,n);
for i=n-1:-1:1
   x(i) = (y(i) - dot(U(i,i+1:end),x(i+1:end)))/U(i,i);
end
end

