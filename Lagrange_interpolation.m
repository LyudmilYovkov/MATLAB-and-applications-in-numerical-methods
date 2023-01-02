%====================================%
% Lagrange interpolation using the   %
% built-in functions polyfit() and   %
% polyval() upon experimental data.  %
%====================================%
% Test function: u(x) = sqrt(x).     %
%====================================%

clear; clc;

%=========%
% Data.   %
%=========%
xdata = [100,121,144];
ydata = [10,11,12];

%=================%
% Interpolation.  %
%=================%
n = length(xdata);
p = polyfit(xdata,ydata,n-1);

syms x
s = 0;
for k = 1 : n
    s = s + p(k) * x^(n-k);
end
s = expand(s);
simplify(s);

display('===========================================')
display('Lagrange interpolating polynomial: ')
pretty(s)
display('===========================================')

%===================================%
% Plot the interpolating polynomial %
% against the experimental data.    %
%===================================%
xi = min(xdata) : 0.01 : max(xdata);
yi = subs(s,x,xi);
yExact = sqrt(xi);
% % %
figure(1)
plot(xi,yi,'b','LineWidth',3)
hold on
grid on
plot(xi,yExact,'g--','LineWidth',3)
plot(xdata,ydata,'ko','LineWidth',2)
axis([min(xi), max(xi), min(double(yi))-0.5, max(double(yi))+0.5])
legend('Polynomial','Exact','Data')

%=========================%
% Estimating the error.   %
%=========================%
r = input('Enter a node (r) for interpolating at: r = ');
u = sqrt(x);
up = diff(u,x);
upp = diff(up,x);
uppp = diff(upp,x);
yr = polyval(p,r); % approximated value
ur = subs(u,x,r); % exact value
upppVals = subs(uppp,x,xi);
maxUPPP = max(upppVals);

prod = 1;
for k = 1 : n
    prod = prod * (r - xdata(k));
end

err = maxUPPP / (factorial(n)) * prod;
err = abs(double(err));
display(['Error: ', num2str(err)])
