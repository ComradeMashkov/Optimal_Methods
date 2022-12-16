clc
clear variables
close all

fnc = @(x) -4 * x(1) * x(2) + 3 * x(1) ^ 2 + 6 * x(2) ^ 2 + 8 * sqrt(5) * x(1) + 4 * sqrt(5) * x(2) + 36;
eps = 0.0000005;
x0 = [-sqrt(5); 0];
k = 0;
n = 2;
Kmax = 10 ^ 6;
path = NaN(3, Kmax);

h = 0.001;
xk = ones(n, 1);
xk1 = x0;
while abs(fnc(xk) - fnc(xk1)) > eps && norm(xk - xk1) > eps && k < Kmax
    k = k + 1;
    path(1, k) = xk1(1);
    path(2, k) = xk1(2);
    path(3, k) = fnc(xk1);
    xk = xk1;  
    grad = -([fnc([xk(1) + h; xk(2)]); fnc([xk(1); xk(2) + h])] - [fnc([xk(1) - h; xk(2)]); fnc([xk(1); xk(2) - h])]) / 2 / h;
    f = @(l) fnc(xk + l * grad);
    lambda =  argmin(f, 0, 2, eps, Kmax);
    xk1 = xk + lambda * grad;
end
Xmin = xk1;
path(1, k + 1) = xk1(1);
path(2, k + 1) = xk1(2);
path(3, k + 1) = fnc(xk1);
Fmin = fnc(Xmin);

fprintf('Number of iterations k = %d \n', k);
fprintf('Points of local min xmin1 = %f, xmin2 = %f \n', Xmin(1), Xmin(2));
fprintf('Function min Fmin = %f \n', Fmin);

hold on
[X, Y] = meshgrid(-10 : 0.2 : 2);
Z = -4 .* X .* Y + 3 .* X .^ 2 + 6 .* Y .^ 2 + 8 .* sqrt(5) .* X + 4 .* sqrt(5) .* Y + 36;
contour(X, Y, Z)
grid on
grid minor
plot3(path(1, :), path(2, :), path(3, :),'r', 'LineWidth', 1)
plot3(path(1, k + 1), path(2, k + 1), path(3, k + 1),'y*', 'LineWidth', 1)
title('Contour')
hold off

function arg = argmin(f, a, b, eps, Kmax)
        fi = (1 + sqrt(5)) / 2;
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