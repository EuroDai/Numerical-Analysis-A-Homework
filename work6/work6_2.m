clear
clc

format long

x_0 = [1.6; 1.2];
tol = 10e-6;

F_0 = F(x_0);
dF_0 = [[2*1.6, 2*1.2];
    [2*1.6, -2*1.2]];
B_0 = inv(dF_0);

x_prev = x_0;
x_next = x_prev - B_0* F_0;
disp(x_next)
n = 1;
disp(n)

F_prev = F_0;
F_next = F(x_next);
p = x_next - x_prev;
q = F_next - F_prev;
B_prev = B_0;
B_next = B_prev + (p - B_prev * q)*(p.')*B_prev/((p.')*B_prev*q);
 
while norm(p) >= tol && norm(F_next) >= tol
    x_prev = x_next;
    x_next = x_prev - B_next*F(x_prev);
    p = x_next - x_prev;
    q = F(x_next) - F(x_prev);
    disp(x_next)
    B_prev = B_next;
    B_next = B_prev + (p - B_prev * q)*(p.')*B_prev/((p.')*B_prev*q);
    n = n+1;
    disp(n)
end


function F_c=F(x)
F_c = [(x(1)^2+x(2)^2-4);(x(1)^2-x(2)^2-1)];
end


