clc
clear all

%% Инициализация

plotCircle = @(center,radius)...
    plot(...
    center(1) + radius*cos(0:0.001:2*pi),...
    center(2) + radius*sin(0:0.001:2*pi),...
    'lineWidth',1.4);

a = -10;
b = 10;
radius = 1;
x0 = [1;0];
xNow = [0;0];
xNext = x0;
eps = 0.001;
approximations = x0';
steps = 0;

%% Двумерная минимизация методом проекции градиента

while norm(xNext-xNow) >= eps
    xNow = xNext;
    
    pNow = -gradient(xNow);
    alphaNow = goldenRatioMinimize(xNow,pNow,a,b,eps);
    xNextStep = xNow + alphaNow * pNow;
    
    if isPointInsideSet(xNextStep,x0,radius)
        xNext = xNextStep;
    else
        approximations = vertcat(approximations,xNextStep');
        xNext = x0 + (xNextStep - x0) * radius / norm(xNextStep - x0);
    end
    
    approximations = vertcat(approximations,xNext');
    steps = steps + 1;
end

%% Представление результатов

% Вывод количества шагов прибилижения
fprintf('Количество шагов приближения: %d\n',steps);
% Вывод точки минимума
fprintf('Проекция точки минимума на множество U: (%f,%f)\n',xNext(1),xNext(2));

plotCircle(x0,radius);
hold on;

% Вывод последовательности приближений в виде графика
plot(...
    approximations(:,1),...
    approximations(:,2),...
    'r','LineWidth',1.25);
title('Последовательность приближений к точке минимума');
xlabel('x1');
ylabel('x2');
grid on;
hold on;

% Вывод линий уровня
x = -5:0.01:1.25;
[X,Y] = meshgrid(x);
Z = -4*X.*Y + 6*X.^2 + 3*Y.^2 ...
    + 4*sqrt(5)*X + 8*sqrt(5)*Y + 22;
contour(X, Y, Z)


%% Целевая двумерная функция

function f = f(x)
f = -4*x(1)*x(2)+6*x(1)^2+3*x(2)^2+...
    4*sqrt(5)*x(1)+8*sqrt(5)*x(2)+22;
end

%% Значение градиента целевой функции в точке

function grad = gradient(x)
grad = zeros(2,1);
grad(1) = -4*x(2)+12*x(1)+4*sqrt(5);
grad(2) = -4*x(1)+6*x(2)+8*sqrt(5);
end

%% Вспомогательная одномерная функция

function F = F(alpha,x,p)
F = f(x+alpha*p);
end

%% Модифицированный под задачу метод золотого сечения

function [alphaMin] = goldenRatioMinimize(x,p,a,b,eps)
alpha1Factor = (3-sqrt(5))/2;
alpha2Factor = (sqrt(5)-1)/2;
doubleEps = 2*eps;

while true
    alpha1 = a + alpha1Factor*(b-a);
    alpha2 = a + alpha2Factor*(b-a);
    
    if F(alpha1,x,p) <= F(alpha2,x,p)
        b = alpha2;
    else
        a = alpha1;
    end
    
    if b - a < doubleEps
        alphaMin = (a+b)/2;
        break;
    end
end

end

%% Проверка принадлежности точки допустимому множеству 

function value = isPointInsideSet(x,center,radius)
value = (x(1) - center(1))^2 + (x(2) - center(2))^2 <= radius ^ 2;
end