clc
clear variables
close all

fnc = @(x) 4 * x(1) * x(2) + 7 * x(1) * x(1) + 4 * x(2) * x(2) + 6 * sqrt(5) * x(1) - 12 * sqrt(5) * x(2) + 51;
eps = 0.05;
x0 = [0; sqrt(5)];
k = 0;
n = 2;
Kmax = 10 ^ 6;
path = NaN(3, Kmax);

h = 0.001;
xk = ones(n, 1);
xk1 = x0;
Hk = [14, 4; 4, 8];
antigrad1 = -([fnc([xk1(1) + h; xk1(2)]); fnc([xk1(1); xk1(2) + h])] - [fnc([xk1(1) - h; xk1(2)]); fnc([xk1(1); xk1(2) - h])]) / 2 / h;

while abs(fnc(xk) - fnc(xk1)) > eps || norm(xk - xk1) > eps && k < Kmax
    k = k + 1;
    path(1, k) = xk1(1);
    path(2, k) = xk1(2);
    path(3, k) = fnc(xk1);
    antigrad = antigrad1;
    xk = xk1;
    xk1 = xk + inv(Hk) * antigrad;
    antigrad1 = -([fnc([xk1(1) + h; xk1(2)]); fnc([xk1(1); xk1(2) + h])] - [fnc([xk1(1) - h; xk1(2)]); fnc([xk1(1); xk1(2) - h])]) / 2 / h;  
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
plot3(path(1, :), path(2, :), path(3, :),'r', 'LineWidth', 1.5);
plot3(path(1, k + 1), path(2, k + 1), path(3, k + 1),'g*', 'LineWidth', 1.5);
hold off
