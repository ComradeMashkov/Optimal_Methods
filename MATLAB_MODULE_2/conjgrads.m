clc
clear variables
close all

fnc = @(x) 4 * x(1) * x(2) + 7 * x(1) * x(1) + 4 * x(2) * x(2) + 6 * sqrt(5) * x(1) - 12 * sqrt(5) * x(2) + 51;

syms x y
f = 4 * x * y + 7 * x * x + 4 * y * y + 6 * sqrt(5) * x - 12 * sqrt(5) * y + 51;
H = hessian(f, [x, y]);
Q = double(H);

eps = 0.05;
x0 = [0; sqrt(5)];
k = 0;
n = 2;
Kmax = 10 ^ 6;
path = NaN(3, Kmax);

xk1 = x0;
bk = 0;
pk = 0;
while norm(antigrad(xk1, fnc)) > eps && k < Kmax
    k = k + 1;
    path(1, k) = xk1(1);
    path(2, k) = xk1(2);
    path(3, k) = fnc(xk1);
    xk = xk1;
    pk = antigrad(xk, fnc) + bk * pk;
    bk = dot(Q * (-antigrad(xk, fnc)), pk) / dot(Q * pk, pk);
    f = @(l) fnc(xk + l * pk);
    lambda =  argmin(f, 0, 2, eps, Kmax);
    xk1 = xk + lambda * pk;
end

Xmin = xk1;
path(1, k + 1) = xk1(1);
path(2, k + 1) = xk1(2);
path(3, k + 1) = fnc(xk1);
Fmin = fnc(Xmin);

fprintf('Число итераций = %d \n', k);
fprintf('Точка минимума = [%e,  %e] \n', Xmin(1), Xmin(2));
fprintf('Значение функции в точке = %e \n', Fmin);

hold on
[X, Y] = meshgrid(-10 : 0.2 : 10);
Z = 4 .* X .* Y + 7 .* X .* X + 4 .* Y .* Y + 6 .* sqrt(5) .* X - 12 .* sqrt(5) .* Y + 51;
contour(X, Y, Z, 'LineWidth', 1.5);
grid on
grid minor
plot3(path(1, :), path(2, :), path(3, :),'r', 'LineWidth', 1.5)
plot3(path(1, k + 1), path(2, k + 1), path(3, k + 1),'g*', 'LineWidth', 1.5)
hold off

function arg = argmin(f, a, b, eps, Kmax)
        phi = (1 + sqrt(5)) / 2;
        kk = 0;
        x1 = b - (b - a) / phi;
        x2 = a + (b - a) / phi;   
        y1 = f(x1);
        y2 = f(x2);
        while abs(b - a) > eps && kk < Kmax
            kk = kk + 1;
            if y1 >= y2
                a = x1;
                x1 = x2;
                x2 = a + (b - a) / phi;
                y1 = y2;
                y2 = f(x2);
            else
                b  = x2;
                x2 = x1;
                x1 = b - (b - a) / phi;
                y2 = y1;
                y1 = f(x1);
            end
        end
        arg = (a + b) / 2;
end
    
function antigrad = antigrad(xk, fnc)
    h = 0.001;
    antigrad = -([fnc([xk(1) + h; xk(2)]); fnc([xk(1); xk(2) + h])] - [fnc([xk(1) - h; xk(2)]); fnc([xk(1); xk(2) - h])]) / 2 / h;
end