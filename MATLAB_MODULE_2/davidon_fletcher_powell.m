clc
clear variables
close all

fnc = @(x) 4 * x(1) * x(2) + 7 * x(1) * x(1) + 4 * x(2) * x(2) + 6 * sqrt(5) * x(1) - 12 * sqrt(5) * x(2) + 51;
eps = 0.005;
x0 = [0; sqrt(5)];
k = 0;
n = 2;
Kmax = 10 ^ 6;
path = NaN(3, Kmax);

h = 0.001;
xk = ones(n, 1);
xk1 = x0;
Hk1 = [14, 4; 4, 8];
antigrad1 = -([fnc([xk1(1) + h; xk1(2)]); fnc([xk1(1); xk1(2) + h])] - [fnc([xk1(1) - h; xk1(2)]); fnc([xk1(1); xk1(2) - h])]) / 2 / h;

while abs(fnc(xk) - fnc(xk1)) > eps && norm(xk - xk1) > eps && k < Kmax
    k = k + 1;
    path(1, k) = xk1(1);
    path(2, k) = xk1(2);
    path(3, k) = fnc(xk1);
    antigrad = antigrad1;
    Hk = Hk1;
    xk = xk1;
    pk = inv(Hk) * antigrad;
    f = @(l) fnc(xk + l * pk);
    lambda = argmin(f, 0, 2, eps, Kmax);
    xk1 = xk + lambda * pk;
    antigrad1 = -([fnc([xk1(1) + h; xk1(2)]); fnc([xk1(1); xk1(2) + h])] - [fnc([xk1(1) - h; xk1(2)]); fnc([xk1(1); xk1(2) - h])]) / 2 / h;
    sk = xk1 - xk;
    yk = -(antigrad1 - antigrad);         
    Hk1 = Hk - (Hk * sk * sk' * Hk') / (sk' * Hk * sk) + (yk * yk') / (yk' * sk);    
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