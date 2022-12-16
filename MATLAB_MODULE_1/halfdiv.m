clc
clear variables
close all


fnc = @(x) x^2 + 3 * x * (log(x) - 1);

a = 0.5;
b = 1;
eps = 0.05;

[Fmin, Xmin, k] = half_division(fnc, a, b, eps);
fprintf('Xmin = %e, Fmin = %e, k = %d\n', Xmin, Fmin, k);

X = 0 : 0.01 : 10;
Y = X.^2 + 3 .* X .* (log(X) - 1);
plot(X, Y, 'LineWidth', 1.3, 'DisplayName', 'f(x)');
hold on;
plot(Xmin, Fmin, 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Fmin');
grid on; grid minor;
xlabel('x'); ylabel('f(x)');
title('Half Division Minimization');
legend show;

function [Fmin, Xmin, k] = half_division(fnc, a, b, eps)
    k = 1;
    while abs(b - a) > eps
        xmin = (a + b) / 2;
        if fnc(xmin - eps) < fnc(xmin + eps)
            b = xmin;
        else
            a = xmin;
        end
        k = k + 1;
    end
    Xmin = (a + b) / 2;
    Fmin = fnc(Xmin);
end