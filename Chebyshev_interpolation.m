%==========================================%
% We use the zeros of the Chebyshev        %
% polynomials to interpolate experimental  %
% data with minimal possible error.        %    
%==========================================%
clear; clc;
%===================================%
% Data.                             % 
%===================================%
a = 0; b = 1;
n = input('Enter n = ');
i = 0 : n;
xc = (a + b) / 2 + (b - a) / 2 * ...
    cos( (2 * i + 1) * pi / (2 * n + 2) );
yc = sqrt(xc);
%=====================================%
% Using polifit() built-in function.  %
%=====================================%
p = polyfit(xc,yc,n);
%===================================%
% Creating finer spatial grid.      %
%===================================%
xx = a : 0.001 : b;
yy = polyval(p,xx);
ff = sqrt(xx);
%============%
% Plot.      %
%============%
figure(1)
% % %
subplot(1,2,1)
% % %
plot(xx,yy,'b','LineWidth',3)
hold on
grid on
plot(xx,ff,'g--','LineWidth',3)
plot(xc,yc,'bo','LineWidth',2)
set(gca,'FontSize',14)
xlabel('x')
ylabel('y')
legend('\bf{Polynomial p(x)}', ...
       '\bf{Exact f(x)}', ...
       '\bf{Data}')
title('\bf{Chebyshev interpolation}')
% % %
subplot(1,2,2)
% % %
plot(xx,abs(yy-ff),'r:','LineWidth',3)
hold on
grid on
set(gca,'FontSize',14)
xlabel('x')
ylabel('y')
legend('\bf{Error = |f(x)-p(x)|}')
title('\bf{Error at Chebyshev interpolation}')
%=======================%
% Approximation error.  %
%=======================%
abs_err = abs(yy - ff);
rel_err = max(abs_err) ./ (max(ff)) .* 100;
display('====================================')
display('Maximal absolute error: ')
display(max(abs_err))
display('Relative error: ')
display(rel_err)
display('====================================')
