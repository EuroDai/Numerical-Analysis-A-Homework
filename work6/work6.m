clc; clear; close all;

%% ====== 原始矩阵 ======
A = [9 4.5 3;
    -56 -28 -18;
     60 30 19];

A1 = A;  A1(3,3) = 18.95;   % a33减小
A2 = A;  A2(3,3) = 19.05;   % a33增大

%% ====== 幂迭代======
function [eigvals] = power_deflation(A)
    n = size(A,1);
    eigvals = zeros(n,1);
    for k = 1:n
        % 初始化随机向量
        x = randn(n,1);
        x = x / norm(x);
        lambda_old = 0;

        % 幂迭代
        for i = 1:100000
            x_new = A * x;
            x_new = x_new / norm(x_new);
            lambda_new = x_new' * A * x_new;
            if abs(lambda_new - lambda_old) < 1e-12
                break;
            end
            lambda_old = lambda_new;
            x = x_new;
        end
        eigvals(k) = lambda_new;

        % 左特征向量估计 (A' 的主特征向量)
        y = randn(n,1);
        y = y / norm(y);
        lambda_old2 = 0;
        for j = 1:10000
            y_new = A' * y;
            y_new = y_new / norm(y_new);
            lambda_new2 = y_new' * A' * y_new;
            if abs(lambda_new2 - lambda_old2) < 1e-12
                break;
            end
            lambda_old2 = lambda_new2;
            y = y_new;
        end

        % deflation 秩一消去：A <- A - λ*v*u'/(u'*v)
        denom = y' * x;
        if abs(denom) < 1e-14, break; end
        A = A - lambda_new * (x * y') / denom;
    end
end

%% ====== 计算特征值 ======
eig_A  = power_deflation(A);
eig_A1 = power_deflation(A1);
eig_A2 = power_deflation(A2);

%% ====== 与MATLAB内置结果比较 ======
eig_ref  = eig(A);
eig_ref1 = eig(A1);
eig_ref2 = eig(A2);

disp('===== 幂迭代计算结果 =====');
disp('原A的特征值:'),  disp(eig_A);
disp('a33=18.95时特征值:'),  disp(eig_A1);
disp('a33=19.05时特征值:'),  disp(eig_A2);

disp('===== MATLAB eig()参考结果 =====');
disp('原A特征值:'),  disp(eig_ref);
disp('a33=18.95:'), disp(eig_ref1);
disp('a33=19.05:'), disp(eig_ref2);

