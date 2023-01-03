%===========================================%
% Approximation of experimental data        %
% using different types of spline curves.   %
%===========================================%

clear; clc;

%=========================================================%
% Generating experimental data over the interval [a; b]   %
% from different test functions:                          %
% f_test_1(x) = 1 ./ (1 + 25 * x.^2);                     %
% f_test_2(x) = 1 ./ (5 + x.^2 + x.^3);                   %
% f_test_3(x) = sqrt(x.^2 + 1) ./ (1 + sin(x).^2) + ...   %
%               1 ./ (x.^2 + 3);                          %
% f_test_4(x) = sin(x.^2) + log(x - sin(x));              %
%=========================================================%
a = 1; b = 3;
xdata = a : 0.1 : b;
% ydata = 1 ./ (1 + 25 * xdata.^2);
% ydata = 1 ./ (5 + xdata.^2 + xdata.^3);
% ydata = sqrt(xdata.^2 + 1) ./ (1 + sin(xdata).^2) + 1 ./ (xdata.^2 + 3);
ydata = sin(xdata.^2) + log(xdata - sin(xdata));

%=================================%
% Introducing a finer grid for    %
% tabulating and a smoother plot. %
%=================================%
xi = a : 0.01 : b;
% yi = 1 ./ (1 + 25 * xi.^2);
% yi = 1 ./ (5 + xi.^2 + xi.^3);
% yi = sqrt(xi.^2 + 1) ./ (1 + sin(xi).^2) + 1 ./ (xi.^2 + 3);
yi = sin(xi.^2) + log(xi - sin(xi));

%=====================================%
% Using interp1() built-in function   %
% with different parameter options.   %
%=====================================%
yLinear = interp1(xdata,ydata,xi,'linear');
ySpline = interp1(xdata,ydata,xi,'spline');
yCubic  = interp1(xdata,ydata,xi,'cubic');

%======================================%
% Comparing the interpolation results  %
% in a graphical way.                  % 
%======================================%
figure(1)
% % %
subplot(2,2,1)
% % %
plot(xi,yLinear,'b','LineWidth',3)
hold on
grid on
plot(xi,yi,'cyan--','LineWidth',3)
plot(xdata,ydata,'ko','LineWidth',1)
legend('option = linear','Exact','Data')
% % %
subplot(2,2,2)
% % %
plot(xi,ySpline,'m','LineWidth',3)
hold on
grid on
plot(xi,yi,'cyan--','LineWidth',3)
plot(xdata,ydata,'ko','LineWidth',1)
legend('option = spline','Exact','Data')
% % %
subplot(2,2,3)
% % %
plot(xi,yCubic,'g','LineWidth',3)
hold on
grid on
plot(xi,yi,'cyan--','LineWidth',3)
plot(xdata,ydata,'ko','LineWidth',1)
legend('option = cubic','Exact','Data')
% % %
subplot(2,2,4)
% % %
plot(xi,yi,'cyan','LineWidth',3)
hold on
grid on
plot(xdata,ydata,'ko','LineWidth',1)
legend('Exact','Data')
  
%===================================%
% Calculating the absolute error    %
% for every approximaton.           %
%===================================%
abs_err_Linear = abs(yi - yLinear);
abs_err_Spline = abs(yi - ySpline);
abs_err_Cubic  = abs(yi - yCubic);

%===================================%
% Calculating the relative error    %
% for every approximation.          %
%===================================%
rel_err_Linear = max(abs_err_Linear) ./ abs(max(yi)) .* 100;
rel_err_Spline = max(abs_err_Spline) ./ abs(max(yi)) .* 100;
rel_err_Cubic  = max(abs_err_Cubic)  ./ abs(max(yi)) .* 100;

%===================================%
% Printing important information.   %
%===================================%
display('============================================')
display(['Maximal absolute error (linear): ', num2str(max(abs_err_Linear))])
display(['Maximal relative error (linear): ', num2str(rel_err_Linear)])
display('============================================')
% % %
display('============================================')
display(['Maximal absolute error (spline): ', num2str(max(abs_err_Spline))])
display(['Maximal relative error (spline): ', num2str(rel_err_Spline)])
display('============================================')
% % %
display('============================================')
display(['Maximal absolute error (cubic): ', num2str(max(abs_err_Cubic))])
display(['Maximal relative error (cubic): ', num2str(rel_err_Cubic)])
display('============================================')
