clc
clear variables
close all


fnc = @(x) x^2 + 3 * x * (log(x) - 1);

a = 0.5;
b = 1;
eps = 0.05;
Kmax = 1e6;

function [Fmin, Xmin, k] = golden_ratio(fnc, a, b, Kmax)
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
        if k > Kmax
            break;
        end
    end
    Xmin = (a + b) / 2;
    Fmin = fnc(Xmin);
end