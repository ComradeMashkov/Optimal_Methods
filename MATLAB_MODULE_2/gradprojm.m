clc
clear variables
close all

% Init

plotCircle = @(center, radius) plot(center(1) + radius*cos(0:0.001:2*pi), center(2) + radius*sin(0:0.001:2*pi), 'lineWidth',1.4);

a = -10;
b = 10;
radius = 1;
x0 = [1;0];
x_current = [0;0];
x_next = x0;
eps = 0.001;
approximations = x0';
steps = 0;

% 2 dim optimization (gradients projection)

while norm(x_next-x_current) >= eps
    x_current = x_next;
    
    pNow = -gradient(x_current);
    alphaNow = goldenRatioMinimize(x_current,pNow,a,b,eps);
    xNextStep = x_current + alphaNow * pNow;
    
    if isPointInsideSet(xNextStep,x0,radius)
        x_next = xNextStep;
    else
        approximations = vertcat(approximations,xNextStep');
        x_next = x0 + (xNextStep - x0) * radius / norm(xNextStep - x0);
    end
    
    approximations = vertcat(approximations,x_next');
    steps = steps + 1;
end

% Output

fprintf('Количество шагов приближения: %d\n', steps);
fprintf('Проекция точки минимума на множество U: (%f,%f)\n', x_next(1), x_next(2));

plotCircle(x0, radius);
hold on;

plot(approximations(:, 1), approximations(:, 2), 'r', 'LineWidth', 1.25);
title('Последовательность приближений к точке минимума');
xlabel('x1');
ylabel('x2');
grid on;
hold on;

% Contour lines
x = -5:0.01:1.25;
[X,Y] = meshgrid(x);
Z = -4*X.*Y + 3*X.^2 + 6*Y.^2 ...
    + 8*sqrt(5)*X + 4*sqrt(5)*Y + 36;
contour(X, Y, Z)


% 2d function

function f = f(x)
f = -4 * x(1) * x(2) + 3 * x(1)^2 + 6 * x(2)^2 + 8 * sqrt(5) * x(1) + 4 * sqrt(5) * x(2) + 36;
end

% Gradient value in point

function grad = gradient(x)
grad = zeros(2, 1);
grad(1) = -4 * x(2) + 6 * x(1) + 8 * sqrt(5);
grad(2) = -4 * x(1) + 12 * x(2) + 4 * sqrt(5);
end

% Additional 1d function

function F = F(alpha, x, p)
F = f(x + alpha * p);
end

% Golden Ration minimization

function [alphaMin] = goldenRatioMinimize(x, p, a, b, eps)
alpha1Factor = (3-sqrt(5))/2;
alpha2Factor = (sqrt(5)-1)/2;
doubleEps = 2*eps;

while true
    alpha1 = a + alpha1Factor*(b - a);
    alpha2 = a + alpha2Factor*(b - a);
    
    if F(alpha1,x,p) <= F(alpha2, x, p)
        b = alpha2;
    else
        a = alpha1;
    end
    
    if b - a < doubleEps
        alphaMin = (a + b)/2;
        break;
    end
end

end

% Belongs to variety check

function value = isPointInsideSet(x, center, radius)
value = (x(1) - center(1))^2 + (x(2) - center(2))^2 <= radius ^ 2;
end