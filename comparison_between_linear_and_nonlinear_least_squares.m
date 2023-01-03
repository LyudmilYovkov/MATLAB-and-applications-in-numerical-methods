%===================================%
% Comparison between linear and     %
% non-linear least squares method   %
% applied to one and same           %
% experimental data set.            %
% In the case of non-linear least   %
% suqares method cftool has been    %
% used.                             %
%===================================%

clear all; clc;

%=======================%
% Experimental data.    %
%=======================%
xdata = [0.05, 0.13, 0.20, 0.27, 0.30, 0.33, 0.38, 0.41, 0.47, ...
        0.50, 0.57, 0.64, 0.70, 0.78, 0.86, 0.91, 1.00];
ydata = [1.02, 0.95, 0.81, 0.70,  0.610, 0.490, 0.47, 0.40, 0.34, ...
        0.28, 0.27, 0.26, 0.258, 0.261, 0.272, 0.22, 0.15];

%========================%
% Applying linear LSM.   %
%========================%
n = 7; % degree of the approximating polynomial
poly = polyfit(xdata,ydata,n);

%==================================%
% Finer grid for smoother plot.    %
%==================================%
xFinerGrid = min(xdata) : 0.01 : max(xdata);

%==========================================%
% Tabulating the approximating polynomial  %
% at the points of the finer grid.         %
%==========================================%
yFinerGrid_Linear_LSM = polyval(poly,xFinerGrid);

syms x y_Poly_LSM
y_Poly_LSM = 0;
for i = 1 : length(poly)
    y_Poly_LSM = ...
        y_Poly_LSM + double(poly(i)) * x^(double(length(poly))-i);
end

display('======================================')
display('Linear least squares')
display('======================================')
yPoly = expand(y_Poly_LSM);
simplify(y_Poly_LSM);
pretty(y_Poly_LSM)
    
%========================================%
% Applying non-linear LSM.               %
% Approximating function:                %
% f(x) = a * exp(- [(x-b)/c].^2 )        %
%========================================%
% a, b, c - coefficient of the           %
% exponential approximating polynomial   %
% with 3 terms. These coefficients are   %
% evaluated by cftool and taken directly %
% from there.                            %
%========================================%
a = [0.1843, 0.1585, 1.0237];
b = [0.8847, 0.6475, 0.0277];
c = [0.1933, 0.2378, 0.3639];

syms y_Nonlin_LSM
y_Nonlin_LSM = 0;
for k = 1 : 3
    y_Nonlin_LSM = y_Nonlin_LSM + a(k) * exp(-((x-b(k))./c(k)).^2);
end

sNonlin = expand(y_Nonlin_LSM);
simplify(y_Nonlin_LSM);

display('======================================')
display('Nonlinear least squares')
display('======================================')
pretty(y_Nonlin_LSM)
    
%==========================================%
% Tabulating the approximating non-linear  %
% exponential function at the points       %
% of the finer grid.                       %
%==========================================%
yFinerGrid_Nonlin_LSM = subs(y_Nonlin_LSM,x,xFinerGrid);

%===============================%
% Graphical comparison between  %
% linear and non-linear LSM.    %
%===============================%
figure(1)
% % %
subplot(2,2,1)
% % %
plot(xFinerGrid,yFinerGrid_Linear_LSM,'g','LineWidth',3)
hold on
grid on
plot(xdata,ydata,'bo','LineWidth',2)
axis([min(xdata),max(xdata),min(ydata)-0.03,max(ydata)+0.03])
legend('Linear least squares','Data')
% % %
subplot(2,2,2)
% % %
plot(xFinerGrid,yFinerGrid_Nonlin_LSM,'r','LineWidth',3)
hold on
grid on
plot(xdata,ydata,'bo','LineWidth',2)
axis([min(xdata),max(xdata),min(ydata)-0.03,max(ydata)+0.03])
legend('Nonlinear least squares','Data')
% % %
subplot(2,2,[3,4])
% % %
plot(xFinerGrid,yFinerGrid_Linear_LSM,'b','LineWidth',3)
hold on
grid on
plot(xFinerGrid,yFinerGrid_Nonlin_LSM,'g--','LineWidth',3)
plot(xdata,ydata,'ko','LineWidth',2)
axis([min(xdata),max(xdata),min(ydata)-0.03,max(ydata)+0.03])
legend('Linear least squares','Nonlinear least squares','Data')

%============================================%
% Difference between them by absolute value. %
%============================================%
abs_err = double(abs(yFinerGrid_Nonlin_LSM - yFinerGrid_Linear_LSM));
display('=========================================================')
display('Difference between the linear and the nonlinear squares: ')
display('=========================================================')
display(max(abs_err))
    
%========================================%
% Conclusion: sometimes using nonlinear  %
% squares is as much effective as using  %
% linear ones.                           %
%========================================%
