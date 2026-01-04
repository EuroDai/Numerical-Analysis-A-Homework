experiment()
function experiment()
    % 病态Hilbert矩阵实验
    clc; clear; close all;

    %% 1. 条件数分析
    fprintf('1. 条件数分析\n');
    n_cond_range = 2:16;
    conds = zeros(size(n_cond_range));
    
    for i = 1:length(n_cond_range)
        n = n_cond_range(i);
        H = hilb(n);
        % 1-范数条件数
        conds(i) = cond(H, 1); 
    end
    
    % 绘制条件数变化图
    figure(1);
    semilogy(n_cond_range, conds, '-bo', 'LineWidth', 1.5);
    title('cond(H)_1随n的变化');
    xlabel('矩阵阶数n');
    ylabel('条件数（对数坐标）');
    grid on;

    %% 2. n=6时的求解与误差计算
    fprintf('2. n=6时的求解误差\n');
    n = 6;
    [H, b, x_true] = setup_problem(n);
    
    % (1) LU分解法
    x_lu = solve_lu_manual(H, b);
    err_lu = norm(x_lu - x_true, inf) / norm(x_true, inf);
    
    % (2) Jacobi迭代法
    [x_jac, iter_jac] = solve_jacobi(H, b);
    err_jac = norm(x_jac - x_true, inf) / norm(x_true, inf);
    
    % (3) SOR迭代法
    omegas = [1.0, 1.25, 1.5];
    err_sor = zeros(length(omegas), 1);
    iters_sor = zeros(length(omegas), 1);
    for i = 1:length(omegas)
        [x_s, k_s] = solve_sor(H, b, omegas(i));
        err_sor(i) = norm(x_s - x_true, inf) / norm(x_true, inf);
        iters_sor(i) = k_s;
    end
    
    % (4) 共轭梯度法
    [x_cg, iter_cg] = solve_cg(H, b);
    err_cg = norm(x_cg - x_true, inf) / norm(x_true, inf);
    
    % 输出结果表格
    fprintf('%-15s %-10s %-15s\n', '方法', '迭代次数', '相对误差');
    fprintf('%-15s %-10s %-.4e\n', 'LU分解', '-', err_lu);
    fprintf('%-15s %-10d %-.4e\n', 'Jacobi', iter_jac, err_jac);
    fprintf('%-15s %-10d %-.4e\n', 'SOR(w=1.0)', iters_sor(1), err_sor(1));
    fprintf('%-15s %-10d %-.4e\n', 'SOR(w=1.25)', iters_sor(2), err_sor(2));
    fprintf('%-15s %-10d %-.4e\n', 'SOR(w=1.5)', iters_sor(3), err_sor(3));
    fprintf('%-15s %-10d %-.4e\n', 'CG', iter_cg, err_cg);
    fprintf('\n');

    %% 3. 误差随n的变化
    fprintf('3. 误差随n的变化\n');
    n_max = 200; 
    err_history = zeros(n_max, 4);
    for k = 2:n_max
        [Hk, bk, xk_true] = setup_problem(k);
        
        % LU
        xk_lu = solve_lu_manual(Hk, bk);
        err_history(k, 1) = norm(xk_lu - xk_true, inf) / norm(xk_true, inf);
        
        % Jacobi
        [xk_jac, ~] = solve_jacobi(Hk, bk);
        err_history(k, 2) = norm(xk_jac - xk_true, inf) / norm(xk_true, inf);
        
        % SOR (w=1.25)
        [xk_sor, ~] = solve_sor(Hk, bk, 1.25);
        err_history(k, 3) = norm(xk_sor - xk_true, inf) / norm(xk_true, inf);
        
        % CG
        [xk_cg, ~] = solve_cg(Hk, bk);
        err_history(k, 4) = norm(xk_cg - xk_true, inf) / norm(xk_true, inf);
        
        % 检查误差是否达到100%
        methods = {'LU', 'Jacobi', 'SOR(1.25)', 'CG'};
        for m = 1:4
            if err_history(k, m) >= 1.0 && err_history(k-1, m) < 1.0
                fprintf('当n = %d时，%s方法的误差达到100%%。\n', ...
                    k, methods{m});
            end
        end
    end
    
    % 绘制误差曲线
    figure(2);
    n_axis = 2:n_max;
    semilogy(n_axis, err_history(2:end, 1), '-o', 'DisplayName', 'LU分解'); hold on;
    semilogy(n_axis, err_history(2:end, 2), '-s', 'DisplayName', 'Jacobi');
    semilogy(n_axis, err_history(2:end, 3), '-^', 'DisplayName', 'SOR(w=1.25)');
    semilogy(n_axis, err_history(2:end, 4), '-x', 'DisplayName', 'CG');
    
    yline(1.0, '--r', '误差 100%', 'LineWidth', 2);
    title('相对误差随n的变化');
    xlabel('矩阵阶数n'); ylabel('相对误差（对数坐标）');
    legend('Location', 'southeast');
    grid on;
end

%% 辅助函数定义 
function [H, b, x_true] = setup_problem(n)
    H = hilb(n);
    x_true = ones(n, 1);
    b = H * x_true;
end

% 1. LU分解
function x = solve_lu_manual(A, b)
    n = length(b);
    L = eye(n); U = zeros(n);
    
    for i = 1:n
        for j = i:n
            U(i,j) = A(i,j) - L(i,1:i-1)*U(1:i-1,j);
        end
        for j = i+1:n
            L(j,i) = (A(j,i) - L(j,1:i-1)*U(1:i-1,i)) / U(i,i);
        end
    end
    
    % 前代解
    y = zeros(n, 1);
    for i = 1:n
        y(i) = b(i) - L(i, 1:i-1) * y(1:i-1);
    end
    
    % 回代解
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i,i);
    end
end

% 2. Jacobi
function [x, k] = solve_jacobi(A, b)
    n = length(b);
    x = zeros(n, 1);
    D = diag(A);
    R = A - diag(D);
    max_iter = 5000;
    tol = 1e-6;
    
    for k = 1:max_iter
        x_new = (b - R*x) ./ D;
        if norm(x_new - x, inf) < tol
            x = x_new; return;
        end
        x = x_new;
    end
end

% 3. SOR
function [x, k] = solve_sor(A, b, w)
    n = length(b);
    x = zeros(n, 1);
    max_iter = 5000;
    tol = 1e-6;
    
    for k = 1:max_iter
        x_old = x;
        for i = 1:n
            % sigma = sum(A(i,j) * x(j)) for j != i
            sigma = A(i, :) * x - A(i, i) * x(i); 
            
            % Gauss-Seidel更新值
            x_gs = (b(i) - sigma) / A(i, i);
            
            % SOR加权更新
            x(i) = (1 - w) * x(i) + w * x_gs;
        end
        
        if norm(x - x_old, inf) < tol
            return;
        end
    end
end

% 4. 共轭梯度法
function [x, k] = solve_cg(A, b)
    n = length(b);
    x = zeros(n, 1);
    r = b - A*x;
    p = r;
    max_iter = 200;
    tol = 1e-6;
    
    for k = 1:max_iter
        if norm(r) < tol, return; end
        
        Ap = A * p;
        alpha = (r' * r) / (p' * Ap);
        x = x + alpha * p;
        r_new = r - alpha * Ap;
        
        beta = (r_new' * r_new) / (r' * r);
        p = r_new + beta * p;
        
        r = r_new;
    end
end