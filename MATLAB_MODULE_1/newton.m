clc
clear variables
close all


fnc = @(x) x^2 + 3 * x * (log(x) - 1);
fnc_1 = @(x) 2 * x + 3 * log(x);
fnc_2 = @(x) 2 + 3 / x;

a = 0.5;
b = 1;
x0 = a + (b - a) * rand(1, 1);
eps = 0.005;

[Fmin, Xmin, k] = newton_optimal(fnc, fnc_1, fnc_2, x0, eps);
fprintf('Xmin = %e, Fmin = %e, k = %d\n', Xmin, Fmin, k);

X = 0 : 0.01 : 10;
Y = X.^2 + 3 .* X .* (log(X) - 1);
plot(X, Y, 'LineWidth', 1.3, 'DisplayName', 'f(x)');
hold on;
plot(Xmin, Fmin, 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Fmin');
grid on; grid minor;
xlabel('x'); ylabel('f(x)');
title('Newton Minimization');
legend show;

function [Fmin, Xmin, k] = newton_optimal(fnc, fnc_1, fnc_2, x0, eps)
    k = 1;
    x_curr = fnc_1(x0) / fnc_2(x0);
    while abs(x_curr) >= eps
        x_curr = fnc_1(x0) / fnc_2(x0);
        x0 = x0 - x_curr;
        k = k + 1;
    end
    Xmin = x0;
    Fmin = fnc(x0);
end