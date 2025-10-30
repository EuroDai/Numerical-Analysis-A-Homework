clear
clc

format long
% 定义方程
f1 = @(x, y) 3*x.^2 - y.^2;         % 方程1
f2 = @(x, y) 3*x.*y.^2 - x.^3 - 1;  % 方程2

% 绘图范围
figure;
fimplicit(f1, [-8 8 -8 8], 'LineWidth', 2); hold on;
fimplicit(f2, [-8 8 -8 8], 'r', 'LineWidth', 2);

% 图形设置
legend('3x^2 - y^2 = 0', '3xy^2 - x^3 - 1 = 0', 'Location', 'best');
xlabel('x');
ylabel('y');
axis equal;
grid on;

x_0 = [1; 1];
x_next = Newton_iter(x_0);
x_prev = x_0;
tol = 10e-6;
while norm(x_next-x_prev) >= tol
x_prev = x_next;
x_next = Newton_iter(x_prev);
disp(x_next)
end

function x_next = Newton_iter(x_prev)
x = x_prev(1);
y = x_prev(2);
dF1_x = 3*x;
dF1_y = -2*y;
dF2_x = 3*(y^2) - 3 *x^2;
dF2_y = 6*x*y;
dF = [[dF1_x, dF1_y];
    [dF2_x, dF2_y]];
F = [(3*(x^2)-y^2); 
    (3*x*(y^2) - x^3 -1)];
x_next = [x; y] - dF \ F;

return
end