clear; clc;

A = [10 7  8  7;
     7  5  6  5;
     8  6 10  9;
     7  5  9 10];

B = [ 5   -1    0    0    0;
     -1   4.5  0.2   0    0;
      0   0.2  1   -0.4   0;
      0    0  -0.4  3     1;
      0    0   0    1     3];

ns = [5 6 10];

disp('使用 eig 计算特征值')
disp('eig(A) =')
disp(eig(A))

disp('eig(B) =')
disp(eig(B))

for k = 1:length(ns)
    n = ns(k);
    C = hilb(n);
    disp(['eig(Hilbert(', num2str(n), ')) ='])
    disp(eig(C))
end

disp('使用基本 QR 迭代计算特征值')
tol = 1e-12;
maxIter = 5000;

lambdaA_qr = basic_qr_eig(A, tol, maxIter);
disp('QR eigenvalues of A =')
disp(lambdaA_qr)

lambdaB_qr = basic_qr_eig(B, tol, maxIter);
disp('QR eigenvalues of B =')
disp(lambdaB_qr)

for k = 1:length(ns)
    n = ns(k);
    C = hilb(n);
    lambdaC_qr = basic_qr_eig(C, tol, maxIter);
    disp(['QR eigenvalues of Hilbert(', num2str(n), ') ='])
    disp(lambdaC_qr)
end

function lambda = basic_qr_eig(M, tol, maxIter)
    if nargin < 2
        tol = 1e-12;
    end
    if nargin < 3
        maxIter = 5000;
    end
    Ak = M;
    for k = 1:maxIter
        [Q,R] = qr(Ak);
        Ak = R*Q;
        off = Ak - diag(diag(Ak));
        if norm(off,'fro') < tol
            break;
        end
    end
    lambda = diag(Ak);
end
