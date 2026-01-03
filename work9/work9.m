f = @(x) 1./(1 + 9*x.^2);

n = 10;
k = 1:(n+1);
xk = cos((2*k - 1)*pi/(2*(n+1)));
yk = f(xk);

c = polyfit(xk, yk, n);

xx = linspace(-1, 1, 400);
ff = f(xx);
L10 = polyval(c, xx);

figure;
plot(xx, ff, 'k-', 'LineWidth', 1.5); hold on;
plot(xx, L10, 'r--', 'LineWidth', 1.5);
plot(xk, yk, 'bo', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
legend('f(x)=1/(1+9x^2)', '\tilde{L}_{10}(x)', 'Chebyshev nodes', 'Location', 'Best');
xlabel('x');
ylabel('y');
grid on;
title('f(x) and 10th-degree Lagrange interpolant on Chebyshev nodes');
