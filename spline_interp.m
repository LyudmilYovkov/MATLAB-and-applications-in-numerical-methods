%===================================
% APPROXIMATION VIA SPLINES
%===================================
clear;
clc;
tic
%===================================
% DATA
%===================================
a = 1; b = 3;
x = a : 0.1 : b;
%y = 1./(1+25*x.^2);
%y = 1./(5+x.^2+x.^3);
%y = sqrt(x.^2+1)./(1+sin(x).^2)+1./(x.^2+3);
y = sin(x.^2)+log(x-sin(x));
%===================================
% FINER GRID FOR INTERPOLATING
%===================================
xi = a : 0.01 : b;
%yi = 1./(1+25*xi.^2);
%yi = 1./(5+xi.^2+xi.^3);
%yi = sqrt(xi.^2+1)./(1+sin(xi).^2)+1./(xi.^2+3);
yi = sin(xi.^2)+log(xi-sin(xi));
%===================================
% A FEW TYPES OF SPLINES
%===================================
y1 = interp1(x,y,xi,'linear');
y2 = interp1(x,y,xi,'spline');
y3 = interp1(x,y,xi,'cubic');
%===================================
% PLOTTING SPLINES IN SUBWINDOWS AND 
% COMPARING WITH THE EXACT FUNCTION 
%===================================
figure(1)
%%%
subplot(2,2,1)
%%%
plot(xi,y1,'b','LineWidth',3)
hold on
grid on
plot(xi,yi,'cyan--','LineWidth',3)
plot(x,y,'bo','LineWidth',2)
legend('Linear','Exact','Data')
%%%
subplot(2,2,2)
%%%
plot(xi,y2,'m','LineWidth',3)
hold on
grid on
plot(xi,yi,'cyan--','LineWidth',3)
plot(x,y,'bo','LineWidth',2)
legend('Spline','Exact','Data')
%%%
subplot(2,2,3)
%%%
plot(xi,y3,'r','LineWidth',3)
hold on
grid on
plot(xi,yi,'cyan--','LineWidth',3)
plot(x,y,'bo','LineWidth',2)
legend('Cubic','Exact','Data')
%%%
subplot(2,2,4)
%%%
plot(xi,yi,'cyan','LineWidth',3)
hold on
grid on
plot(x,y,'bo','LineWidth',2)
legend('Exact','Data')
%===================================
% ABSOLUTE ERRORS
%===================================
abs_err1 = abs(yi-y1);
abs_err2 = abs(yi-y2);
abs_err3 = abs(yi-y3);
%===================================
% RELATIVE ERRORS
%===================================
rel_err1 = max(abs_err1)./abs(max(yi)).*100;
rel_err2 = max(abs_err2)./abs(max(yi)).*100;
rel_err3 = max(abs_err3)./abs(max(yi)).*100;
%===================================
% DISPLAY INFORMATION
%===================================
display('============================================')
display('Maximal absolute error (linear): ')
display(max(abs_err1))
display('Maximal relative error (linear): ')
display(rel_err1)
display('============================================')
%%%
display('============================================')
display('Maximal absolute error (spline): ')
display(max(abs_err2))
display('Maximal relative error (spline): ')
display(rel_err2)
display('============================================')
%%%
display('============================================')
display('Maximal absolute error (cubic): ')
display(max(abs_err3))
display('Maximal relative error (cubic): ')
display(rel_err3)
display('============================================')
%===================================
% INTERPOLATION AT INNER POINTS
%===================================
display('============================================')
x0 = input('Enter a number x0 in [a;b] for interpolating: x0 = ');
index = find(abs(xi-x0)<1.0e-03);
display('Linear: ')
display(y1(index(1)))
display('Spline: ')
display(y2(index(1)))
display('Cubic: ')
display(y3(index(1)))
display('============================================')
toc