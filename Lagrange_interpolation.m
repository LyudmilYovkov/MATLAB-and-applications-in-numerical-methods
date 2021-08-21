%===================================
% LAGRANGE INTERPOLATION; USING AN
% EMBEDDED FUNCTIONS POLYFIT() AND
% POLYVAL()
%===================================
clear; clc;
%===================================
% DATA
%===================================
x0 = [100,121,144];
y0 = [10,11,12];
%===================================
% INTERPOLATION
%===================================
n = length(x0);
p = polyfit(x0,y0,n-1);
syms x
s = 0;
for k = 1 : n
    s = s + p(k)*x^(n-k);
end
s = expand(s);
simplify(s);
display('===========================================')
display('Lagrange interpolating polynomial: ')
pretty(s)
display('===========================================')
%===================================
% PLOT THE DATA AND THE POLYNOMIAL
%===================================
xi = min(x0) : 0.01 : max(x0);
yi = subs(s,x,xi);
yExact = sqrt(xi);
%%%
figure(1)
plot(xi,yi,'b','LineWidth',3)
hold on
grid on
plot(xi,yExact,'g--','LineWidth',3)
plot(x0,y0,'bo','LineWidth',2)
axis([min(xi),max(xi),min(double(yi))-0.5,max(double(yi))+0.5])
legend('Polynomial','Exact','Data')
%===================================
% ERROR
%===================================
display('===========================================')
r = input('Enter a node for interpolating: r = ');
u = sqrt(x);
up = diff(u,x);
upp = diff(up,x);
uppp = diff(upp,x);
yr = polyval(p,r); % approximated
ur = subs(u,x,r); % exact
upppVals = subs(uppp,x,xi);
maxUPPP = max(upppVals);
prod = 1;
for k = 1 : n
    prod = prod * (r-x0(k));
end
err = maxUPPP/(factorial(n)) * prod;
err = abs(double(err));
display('Error: ')
display(err)
display('===========================================')