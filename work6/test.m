clear
clc
A = [[4, 2, 1];
    [2, 1, 1];
    [1, 1, -1]];
A_0 = A-eye(3);
B = inv(A);
v_prev = [1; 1; 1];
lambda_prev = lam(v_prev, A);
v_next = rp(v_prev, B);
disp(v_next)
lambda_next = lam(v_next, A);
n=1;
disp(n);
while abs(lambda_prev-lambda_next) >10^(-3)
    v_prev = v_next;
    lambda_prev = lambda_next;
    v_next = rp(v_prev, B);
    disp(v_next)
    lambda_next = lam(v_next, A);
    n = n + 1;
    disp(n);
end
lambda = 

function lambda = lam(v, A)
lambda = transpose(v) * A * v/(transpose(v) * v);
end

function v_next = rp(v, B)
z = B * v;
v_next = z/norm(z);
end