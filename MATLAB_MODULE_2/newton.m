clc
clear variables
close all

fnc = @(x) -4 * x(1) * x(2) + 3 * x(1) ^ 2 + 6 * x(2) ^ 2 + 8 * sqrt(5) * x(1) + 4 * sqrt(5) * x(2) + 36;
eps = 0.05;
x0 = [-sqrt(5); 0];
k = 0;
n = 2;
Kmax = 10 ^ 6;
path = NaN(3, Kmax);

h = 0.001;
xk = ones(n, 1);
xk1 = x0;
Hk = [6, -4; -4, 12];
antigrad1 = -([fnc([xk1(1) + h; xk1(2)]); fnc([xk1(1); xk1(2) + h])] - [fnc([xk1(1) - h; xk1(2)]); fnc([xk1(1); xk1(2) - h])]) / 2 / h;

while abs(fnc(xk) - fnc(xk1)) > eps && norm(xk - xk1) > eps && k < Kmax
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

fprintf('Number of iterations = %d \n', k);
fprintf('Minimum point = [%.3f,  %.3f] \n', Xmin(1), Xmin(2));
fprintf('Function minimum = %.f \n', Fmin);

hold on
[X, Y] = meshgrid(-10 : 0.2 : 2);
Z = -4 .* X .* Y + 3 .* X .^ 2 + 6 .* Y .^ 2 + 8 .* sqrt(5) .* X + 4 .* sqrt(5) .* Y + 36;
contour(X, Y, Z, 'LineWidth', 1.5);
grid on
grid minor
plot3(path(1, :), path(2, :), path(3, :),'r', 'LineWidth', 1.5);
plot3(path(1, k + 1), path(2, k + 1), path(3, k + 1),'g*', 'LineWidth', 1.5);
hold off
