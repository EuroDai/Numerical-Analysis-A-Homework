clear all
clc
format long g

% 数据
A1=[3.01 6.03 1.99; 1.27 4.16 -1.23; 0.987 -4.81 9.34];
A2=[3.00 6.03 1.99; 1.27 4.16 -1.23; 0.990 -4.81 9.34];
b =[1;1;1];

% 计算
x1=cp(A1,b);  % 求解(1)
x2=cp(A2,b);  % 求解(2)
d1=det(A1);   % A1行列式
d2=det(A2);   % A2行列式
k1=cond(A1);  % (1)的条件数

% A1x=b
disp('A1 ='); disp(A1)
disp('b =');  disp(b)
fprintf('det(A1) = %.12g\n', d1)
disp('x1 =');   disp(x1)
fprintf('cond(A1) = %.12g（2-范数条件数）\n\n', k1)

% A2x=b
disp('A2 ='); disp(A2)
disp('b =');  disp(b)
fprintf('det(A2) = %.12g\n', d2)
disp('x2 =');   disp(x2)

% ----- 对比：仅 a11 与 a31 有微小改动 -----
fprintf('\na11 与 a31 有微小差别：\n')
fprintf('Δa11 = %.12g, Δa31 = %.12g\n', A2(1,1)-A1(1,1), A2(3,1)-A1(3,1))
fprintf('||A||∞的相对变化：%.12g\n', norm(A2-A1,inf)/norm(A1,inf))
fprintf('||x||∞的相对变化：%.12g\n', norm(x2-x1,inf)/norm(x1,inf))
fprintf('结论：若 cond(A1) 很大，则极小的系数扰动也会导致解发生显著变化。\n')

% 列主元素消去法
function x=cp(A,b)
n=size(A,1); A=[A b];
for k=1:n-1
    [~,p]=max(abs(A(k:n,k)));
    p=p+k-1;
    A([k p],:)=A([p k],:);
    for i=k+1:n
        A(i,k:n+1)=A(i,k:n+1)-A(i,k)/A(k,k)*A(k,k:n+1);
    end
end
x=zeros(n,1);
for i=n:-1:1
    x(i)=(A(i,end)-A(i,i+1:n)*x(i+1:n))/A(i,i);
end
end
