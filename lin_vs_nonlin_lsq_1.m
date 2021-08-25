%===================================
% COMPARISON BETWEEN LINEAR AND 
% NONLINEAR LEAST SQUARES; USING
% CFTOOL FOR THE NONLINEAR SQUARES
%===================================
clear all
clc
tic
%===================================
% DATA
%===================================
xdat = [0.05, 0.13, 0.2, 0.27, 0.3, 0.33, 0.38, 0.41, 0.47, ...
    0.5, 0.57, 0.64, 0.7, 0.78, 0.86, 0.91, 1];
ydat = [1.02, 0.95, 0.81, 0.7, 0.61, 0.49, 0.47, 0.4, 0.34, ...
    0.28, 0.27, 0.26, 0.258, 0.261, 0.272, 0.22, 0.15];
%===================================
% LINEAR LEAST SQUARES
%===================================
n = 7; 
p = polyfit(xdat,ydat,n);
xi = min(xdat) : 0.01 : max(xdat);
yi = polyval(p,xi);
syms x y
y = 0;
for i = 1 : length(p)
    y = y + p(i) * x^(length(p)-i);
end
display('======================================')
display('Linear least squares')
display('======================================')
y = expand(y);
simplify(y);
pretty(y)
%===================================
% NONLINEAR LEAST SQUARES
%===================================
%===================================
% a, b, c - coefficient of the
% exponential approximating polynomial
% with 3 terms; these coefficients are 
% evaluated by cftool and taken directly 
% from there
%===================================
a = [0.1843, 0.1585, 1.0237];
b = [0.8847, 0.6475, 0.0277];
c = [0.1933, 0.2378, 0.3639];
syms s
s = 0;
for k = 1 : 3
    s = s + a(k) * exp(-((x-b(k))./c(k)).^2);
end
s = expand(s);
simplify(s);
display('======================================')
display('Nonlinear least squares')
display('======================================')
pretty(s)
xx = min(xdat) : 0.01 : max(xdat);
yy = subs(s,x,xx);
%===================================
% PLOT LINEAR VERSUS NONLINEAR
% IN SEPARATE SUBWINDOWS
%===================================
figure(1)
%%%
subplot(2,2,1)
%%%
plot(xi,yi,'g','LineWidth',3)
hold on
grid on
plot(xdat,ydat,'bo','LineWidth',2)
axis([min(xdat),max(xdat),min(ydat)-0.03,max(ydat)+0.03])
legend('Linear least squares','Data')
%%%
subplot(2,2,2)
%%%
plot(xx,yy,'r','LineWidth',3)
hold on
grid on
plot(xdat,ydat,'bo','LineWidth',2)
axis([min(xdat),max(xdat),min(ydat)-0.03,max(ydat)+0.03])
legend('Nonlinear least squares','Data')
%%%
subplot(2,2,[3,4])
%%%
plot(xi,yi,'b','LineWidth',3)
hold on
grid on
plot(xx,yy,'g--','LineWidth',3)
plot(xdat,ydat,'ko','LineWidth',2)
axis([min(xdat),max(xdat),min(ydat)-0.03,max(ydat)+0.03])
legend('Linear least squares','Nonlinear least squares','Data')
%===================================
% Difference between them
%===================================
Diff = double(abs(yy-yi));
display('=========================================================')
display('Difference between the linear and the nonlinear squares: ')
display('=========================================================')
display(max(Diff))
%===================================
% Conclusion: sometimes using nonlinear
% squares is as much effective as using
% linear ones
%===================================
toc