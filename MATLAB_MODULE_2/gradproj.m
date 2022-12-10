clc
clear variables
close all

tic
% Init

plotCircle = @(center, radius) plot(center(1) + radius * cos(0 : 0.001 : 2 * pi), center(2) + radius * sin(0 : 0.001 : 2 * pi), 'lineWidth', 1.4);

a = -10;
b = 10;
radius = 1;
x0 = [1; 0];
x_current = [0; 0];
x_next = x0;
eps = 0.001;
approx = x0';
steps = 0;

% 2 dim optimization (gradients projection)

while norm(x_next-x_current) >= eps
    x_current = x_next;
    
    p_current = -gradient(x_current);
    alpha_current = goldenRatioMinimize(x_current, p_current, a, b, eps);
    x_next_step = x_current + alpha_current * p_current;
    
    if isPointInsideSet(x_next_step, x0, radius)
        x_next = x_next_step;
    else
        approx = vertcat(approx, x_next_step');
        x_next = x0 + (x_next_step - x0) * radius / norm(x_next_step - x0);
    end
    
    approx = vertcat(approx, x_next');
    steps = steps + 1;
end

% Output

fprintf('Number of approximation steps: %d\n', steps);
fprintf('Projection of the minimum point onto the set U: (%e, %e)\n', x_next(1), x_next(2));

plotCircle(x0, radius);
hold on;

plot(approx(:, 1), approx(:, 2), 'r--', 'LineWidth', 1.2);
title('The sequence of approximations to the minimum point', 'FontSize', 14);
xlabel('x_1', 'FontSize', 14);
ylabel('x_2', 'FontSize', 14);
grid on; grid minor;
hold on;

% Contour lines
values = -5 : 0.01 : 2.5;
[X, Y] = meshgrid(values);
Z = -4*X.*Y + 3*X.^2 + 6*Y.^2 + 8*sqrt(5)*X + 4*sqrt(5)*Y + 36;
contour(X, Y, Z, 'LineWidth', 1.3);

toc

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

% Belongs to set check
function value = isPointInsideSet(x, center, radius)
    value = (x(1) - center(1))^2 + (x(2) - center(2))^2 <= radius ^ 2;
end