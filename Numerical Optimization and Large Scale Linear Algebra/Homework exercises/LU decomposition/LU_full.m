clc,clear;

tic
n = 2048;
normi1 = zeros(10,1);
normi2 = zeros(10,1);
for i = 1:10
    %A = eye(n);
    % idx = tril(true(size(A)), -1);
    % A(idx) = -1;
    % A(:,n)=1;
    A = -100 + 200*rand(n,n);
    x_known = ones(n,1);
    %x_known = rand(n,1);
    b = A*x_known;
    [p, q, L, U] = fullLU(A);
    x = solvefullLU(p,q,L,U,b);

    normi1(i) = norm(x-x_known, inf);
    normi2(i) = norm(A*x-b,inf);
end

meanError = mean(normi1);
meanCorrection = mean(normi2);
toc
disp(['The average norm of x error is ',num2str(meanError),' and the average norm of Ax-b correction is ',num2str(meanCorrection)]);

function [p, q, L, U] = fullLU(A)
% fullLU  returns L, U
% matrices and the p and q vectors after 
% full pivoting decomposition of A
LU = A;
n = length(A);
L = eye(n);
p = 1:n;
q = 1:n;
tic
for i = 1:n-1
      LU_tmp = LU(i:end,i:end);
     [~,pos] = max(abs(LU_tmp(:)));
     [p_row, p_col] = ind2sub(size(LU_tmp),pos);
     p([i, p_row+i-1]) = p([p_row+i-1, i]);
     q([i, p_col+i-1]) = q([p_col+i-1, i]);
     LU([p_row+i-1, i], :) = LU([i, p_row+i-1], :);
     LU(:, [p_col+i-1, i]) = LU(:, [i, p_col+i-1]);
     k = i+1:n;
     l = LU(k,i)/LU(i,i);
     LU(k,i:n) =  LU(k,i:n)- l.*LU(i,i:n);
     LU(i+1:end,i) = l;
end
U = triu(LU);
L = tril(LU,-1);
L(1:n+1:end) = ones(n,1);
end

function x = solvefullLU(p,q,L,U,b)
% SolvefullLU returns the
% solution x of the system
% P^{T}LUQ^{T}x = b
n=length(b);
b_LU = b(p);
x = zeros(n,1);
y = zeros(n,1);
z = zeros(n,1);
% front substitution
for i=1:n
    z(i) = b_LU(i) - dot(L(i,:),z);
end
% back substitution
y(n) = z(n)/U(n,n);
 for i=n-1:-1:1
     y(i) = (z(i) - dot(U(i,i+1:end),y(i+1:end)))/U(i,i);
 end
 toc
I = eye(n,n);
I = I(:,q);
x = I*y;
end
