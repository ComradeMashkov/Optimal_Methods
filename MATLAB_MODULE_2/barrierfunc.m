clc
clear variables
close all

tic

global k, k = 1;
global a; 
global r;
a = [1; 0]; r = 1;

fnc = @(x) -4*x(1)*x(2) + 3*x(1)*x(1) + 6*x(2)*x(2) + 8*sqrt(5)*x(1) + 4*sqrt(5)*x(2) + 36;
Fnc = @(x) fnc(x) - 1/k/g(x);

n = 2;
xk1 = zeros(2, 1); 
xk = a; 
eps = 0.0001;
kmax = 1000;
path = NaN(kmax, n + 1);
path(1, :) = [xk', fnc(xk)];
cond = true; x = zeros(kmax, n); x(1, :) = xk';

eps = 0.00001;

while cond
    k = k + 1;
    xk1 = xk;
    
    xk_ = xk;
    xk1_ = xk_ + 2*eps;
    while norm(xk_ - xk1_) > eps  
        xk1_ = xk_;
        for i = 1 : n
            e = zeros(n, 1);
            e(i) = 1;
            f = @(lambda) Fnc(xk_ + lambda*e);
            lambda =  argmin(f, -r, r, eps, kmax);
            xk_ = xk_ + lambda*e;
        end
    end

    xk = xk_;
    path(k, :) = [xk', fnc(xk)];
    x(k, :) = xk';
    if mod(k, 2) == 0
        cond = (norm(x(k, :) - x(k/2, :)) > eps);
    end

end
xmin = xk; fmin = fnc(xmin);
   
figure(1)


fprintf('Number of approximation steps = %d\n', k);
fprintf('xmin = (%e, %e)\n', xmin(1), xmin(2))
fprintf('fmin = %e \n', fmin);

hold on
[X, Y] = meshgrid(-5 : 0.01 : 5);
Z = -4 .* X .* Y + 3 .* X .^ 2 + 6 .* Y .^ 2 + 8 .* sqrt(5) .* X + 4 .* sqrt(5) .* Y + 36;
contour(X, Y, Z)
grid on
grid minor

plot3(path(:, 1), path(:, 2), path(:, 3), 'r--', 'LineWidth', 1.5)
plot3(path(k, 1), path(k, 2), path(k, 3),'b*', 'LineWidth', 3)
xmin_actual = [-2.233000e+00,  4.464583e+00];
plot3(xmin_actual(1), xmin_actual(2), fnc(xmin_actual), 'g*', 'LineWidth', 3)

[XS, YS, ZS] = sphere;
XU = XS * r;
YU = YS * r;
ZU = ZS * r;
surf(XU + a(1), YU + a(2), ZU - 50)

toc

function g = g(x)
    global a; global r;
    g = (x(1) - a(1))^2 + (x(2) - a(2))^2 - r^2;
    if g > 0 
        g = -1e-9;
    end
end

function arg = argmin(f, a, b, eps, kmax)
    t = (sqrt(5) - 1) / 2;
    x1 = a + (1 - t)*(b - a);
    x2 = a + t*(b - a);
    y1 = f(x1);
    y2 = f(x2);
    k = 0;
    while b - a > 2*eps && k < kmax
        k = k + 1;
        if y1 <= y2
            b = x2;
            x2 = x1;
            x1 = a + (1 - t)*(b - a);
            y2 = y1;
            y1 = f(x1);
        else
            a = x1;
            x1 = x2;
            x2 = a + t*(b - a);
            y1 = y2;
            y2 = f(x2);
        end
    end
    arg = (a + b) / 2;
end