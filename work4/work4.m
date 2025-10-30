clear; clc
nList = [6 8 10 12];
tol = 1e-4; maxit = 5e5;
omegas = [1.00 1.25 1.50];

fprintf('%6s %8s %8s %14s\n','n','method','omega','||r||_inf');

for n = nList
    H = hilb(n);
    x_true = ones(n,1); 
    b = H*x_true;
    x0 = zeros(n,1);

    for w = omegas
        [xS,rS] = sor(H,b,w,x0,tol,maxit);
        fprintf('%6d %8s %8.2f %14.3e\n', n, 'SOR', w, rS);
    end
    fprintf('\n');
end

function [x,res] = sor(A,b,omega,x,tol,maxit)
n = length(b);
for k = 1:maxit
    x_old = x;
    for i = 1:n
        s1 = A(i,1:i-1)*x(1:i-1);
        s2 = A(i,i+1:n)*x_old(i+1:n);
        x(i) = (1-omega)*x_old(i) + omega*(b(i)-s1-s2)/A(i,i);
    end
    r = A*x - b;
    res = norm(r,inf);
    if res < tol
        break
    end
end
end