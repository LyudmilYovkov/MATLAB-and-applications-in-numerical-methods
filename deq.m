clear; clc;
fun = @(t,u) [u(2); u(3); sin(t)-u(1)];
[x, y] = ode45(fun, [0, 5], [0; 1; 1]);
hold on 
grid on
axis([0, 5, 0, 6])
set(gca, 'FontName', 'Times', 'FontSize', 14)
title('$$ y'''''' + y = \sin x, \quad y(0)=0, \quad y''(0)=1, \quad y''''(0)=1 $$', ...
    'interpreter', 'latex')
xlabel('$$ x \in [-3; \, 5] $$', 'interpreter', 'latex')
ylabel('$$ y \in [0; \, 6] $$', 'interpreter', 'latex')
y1 = dsolve('D3y + y = sin(x), y(0)=0, Dy(0)=1, D2y(0)=1', 'x');
y1 = sym(y1); % символно, за да можем да използваме simplify() и pretty()
y1 = simplify(y1);
disp('Точното решение е: y(x)=')
pretty(y1)
xx = 0 : 0.05 : 5;
yy = subs(y1, 'x', xx);
plot(x, y(:, 1), 'b', 'LineWidth', 3)
plot(xx, yy, 'c:', 'LineWidth', 3)
legend('\it{Приближеното решение}', '\it{Точното решение}')