clear; clc; close all;

% 原函数
f = @(x) 1./(1 + 9*x.^2);

% 插值节点（改成列向量）
x = (-1:0.2:1)';
y = f(x);
n = length(x) - 1;

%% Newton 插值多项式

dd = y;
for k = 2:n+1
    dd(k:n+1) = (dd(k:n+1) - dd(k-1:n)) ./ (x(k:n+1) - x(1:n-k+2));
end

xx = linspace(-1, 1, 400);
N10_val = dd(end) * ones(size(xx));
for k = n:-1:1
    N10_val = dd(k) + (xx - x(k)) .* N10_val;
end

%% 三次样条
m = n + 1;
h = x(2) - x(1);

A = zeros(m);
rhs = zeros(m,1);

A(1,1) = 1;
A(m,m) = 1;

for i = 2:m-1
    A(i,i-1) = h;
    A(i,i)   = 2*(h + h);
    A(i,i+1) = h;
    rhs(i) = 6 * ( (y(i+1)-y(i))/h - (y(i)-y(i-1))/h );
end

M = A \ rhs;

S3_val = zeros(size(xx));
for j = 1:length(xx)
    xq = xx(j);

    if xq == x(end)
        i = m - 1;
    else
        i = find(xq >= x(1:end-1) & xq < x(2:end), 1, 'last');
    end

    hi = x(i+1) - x(i);
    A1 = x(i+1) - xq;
    B1 = xq - x(i);

    S3_val(j) = M(i)  * A1^3 / (6*hi) + ...
                M(i+1)* B1^3 / (6*hi) + ...
                (y(i)   - M(i)*hi^2/6)   * (A1/hi) + ...
                (y(i+1) - M(i+1)*hi^2/6) * (B1/hi);
end

figure;
plot(xx, f(xx), 'k-', 'LineWidth', 1.5); hold on;
plot(xx, N10_val, 'r--', 'LineWidth', 1.5);
plot(xx, S3_val, 'b-.', 'LineWidth', 1.5);
plot(x, y, 'ko', 'MarkerFaceColor', 'k');

legend('f(x)','Newton N_{10}(x)','Natural spline S_3(x)','nodes',...
       'Location','Best');
xlabel('x');
ylabel('y');
title('f(x) 与 Newton 插值 N_{10}, 自然三次样条 S_3 比较');
grid on;
