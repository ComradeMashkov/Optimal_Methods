clc
clear all
close all
%Метод сопряженных градиентов

fnc = @(x) -4 * x(1) * x(2) + 3 * x(1) ^ 2 + 6 * x(2) ^ 2 + 8 * sqrt(5) * x(1) + 4 * sqrt(5) * x(2) + 36;
Q = [3, -4; -4, 6];
eps = 0.05;
x0 = [-sqrt(5); 0];
k = 0;
n = 2;
Kmax = 10 ^ 6;
path = NaN(3, Kmax);

xk1 = x0;
bk = 0;
pk = 0;
while norm(grad(xk1, fnc)) > eps && k < Kmax
    k = k + 1;
    path(1, k) = xk1(1);
    path(2, k) = xk1(2);
    path(3, k) = fnc(xk1);
    xk = xk1;
    pk = grad(xk, fnc) + bk * pk;
    bk = dot(Q * (-grad(xk, fnc)), pk) / dot(Q * pk, pk);
    f = @(l) fnc(xk + l * pk);
    lambda =  argmin(f, 0, 2, eps, Kmax);
    xk1 = xk + lambda * pk;
end
Xmin = xk1;
path(1, k + 1) = xk1(1);
path(2, k + 1) = xk1(2);
path(3, k + 1) = fnc(xk1);
Fmin = fnc(Xmin);

fprintf('Количество итераций = %d \n', k);
fprintf('Точка минимума = %e,  %e \n', Xmin(1), Xmin(2));
fprintf('Минимум функции = %e \n', Fmin);

hold on
[X, Y] = meshgrid(-10 : 0.2 : 2);
Z = -4 .* X .* Y + 3 .* X .^ 2 + 6 .* Y .^ 2 + 8 .* sqrt(5) .* X + 4 .* sqrt(5) .* Y + 36;
surf(X, Y, Z)
grid on
grid minor
plot3(path(1, :), path(2, :), path(3, :),'r', 'LineWidth', 1)
plot3(path(1, k + 1), path(2, k + 1), path(3, k + 1),'y*', 'LineWidth', 1)
hold off

function arg = argmin(f, a, b, eps, Kmax)
        fi = 1.62;
        kk = 0;
        x1 = b - (b - a) / fi;
        x2 = a + (b - a) / fi;   
        y1 = f(x1);
        y2 = f(x2);
        while abs(b - a) > eps && kk < Kmax
            kk = kk + 1;
            if y1 >= y2
                a = x1;
                x1 = x2;
                x2 = a + (b - a) / fi;
                y1 = y2;
                y2 = f(x2);
            else
                b  = x2;
                x2 = x1;
                x1 = b - (b - a) / fi;
                y2 = y1;
                y1 = f(x1);
            end
        end
        arg = (a + b) / 2;
end
    
function grad = grad(xk, fnc)
    h = 0.001;
    grad = -([fnc([xk(1) + h; xk(2)]); fnc([xk(1); xk(2) + h])] - [fnc([xk(1) - h; xk(2)]); fnc([xk(1); xk(2) - h])]) / 2 / h;
end