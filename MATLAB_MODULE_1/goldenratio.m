clc
clear variables
close all


fnc = @(x) x^2 + 3 * x * (log(x) - 1);

a = 0.5;
b = 1;
eps = 0.05;

[Fmin, Xmin, k] = golden_ratio(fnc, a, b, eps);
fprintf('Xmin = %e, Fmin = %e, k = %d\n', Xmin, Fmin, k);

X = 0 : 0.01 : 10;
Y = X.^2 + 3 .* X .* (log(X) - 1);
plot(X, Y, 'LineWidth', 1.3, 'DisplayName', 'f(x)');
hold on;
plot(Xmin, Fmin, 'o', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Fmin');
grid on; grid minor;
xlabel('x'); ylabel('f(x)');
title('Golden Ratio Minimization');
legend show;

function [Fmin, Xmin, k] = golden_ratio(fnc, a, b, eps)
    k = 1;
    phi = (1 + sqrt(5)) / 2;
    x1 = b - (b - a) / phi;
    x2 = a + (b - a) / phi;
    while abs(b - a) > eps
        y1 = fnc(x1);
        y2 = fnc(x2);
        if y1 >= y2
            a = x1;
            x1 = x2;
            x2 = a + (b - a) / phi;
        else
            b = x2;
            x2 = x1;
            x1 = b - (b - a) / phi;
        end
        k = k + 1;
    end
    Xmin = (a + b) / 2;
    Fmin = fnc(Xmin);
end