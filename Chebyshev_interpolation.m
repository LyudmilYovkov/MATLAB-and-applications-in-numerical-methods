%===================================
% USING THE ZEROF OF THE CHEBISHOV'S
% POLYNOMIALS FOR INTERPOLATING DATA;
% IN THIS AND ONLY IN THIS CASE THE 
% ERROR OF THE APPROXIMATION IS 
% MINIMAL
%===================================
clear;  clc;
%===================================
% DATA
%===================================
a = 0; b = 1;
n = input('Enter n = '); % n=12
i = 0 : n;
xc = (a+b)/2 + (b-a)/2 * ...
    cos((2*i+1)*pi/(2*n+2));
yc = sqrt(xc);
%===================================
% CALLING POLYFIT()
%===================================
p = polyfit(xc,yc,n);
%===================================
% FINER GRID
%===================================
xx = a : 0.001 : b;
yy = polyval(p,xx);
ff = sqrt(xx);
%===================================
% PLOT
%===================================
figure(1)
%%%
subplot(1,2,1)
%%%
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
%%%
subplot(1,2,2)
%%%
plot(xx,abs(yy-ff),'r:','LineWidth',3)
hold on
grid on
set(gca,'FontSize',14)
xlabel('x')
ylabel('y')
legend('\bf{Error = |f(x)-p(x)|}')
title('\bf{Error at Chebyshev interpolation}')
%===================================
% ERROR
%===================================
abs_err = abs(yy-ff);
rel_err = max(abs_err)./(max(ff)).*100;
display('====================================')
display('Maximal absolute error: ')
display(max(abs_err))
display('Relative error: ')
display(rel_err)
display('====================================')