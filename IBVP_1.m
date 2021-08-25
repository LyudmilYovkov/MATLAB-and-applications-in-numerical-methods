tic
%===================================
% BVP for ODE of 2nd order
% u'' - u = x
% u(0) = 0
% u'(1) = (e^2-e+1) / (e^2);
% u_{exact} = e^(-x) * (-1 + e^(2*x) - x*e^(x));
%===================================
clear all
clc
%===================================
% Grid
%===================================
a = 0; b = 1; h = 0.01;
x = a : h : b;
N = length(x);
y = zeros(1,N);
%===================================
% Exact solution
%===================================
ur = exp(-x) .* (-1 + exp(2*x) - x.*exp(x));
%===================================
% Initial condition
%===================================
y(1) = 0;
%===================================
% adiag
%===================================
adiag = zeros(1,N);
adiag(1) = 0;
adiag(end) = -1/h;
for ii = 2 : N-1
    adiag(ii) = 1;
end
%===================================
% bdiag
%===================================
bdiag = zeros(1,N);
bdiag(1) = 1;
bdiag(end) = 1/h + h/2;
for ii = 2 : N-1
    bdiag(ii) = -(h^2+2);
end
%===================================
% cdiag
%===================================
cdiag = zeros(1,N);
cdiag(1) = 0;
cdiag(end) = 0;
for ii = 2 : N-1
    cdiag(ii) = 1;
end
%===================================
% ddiag
%===================================
alpha = (exp(1)^2 - exp(1) + 1) / (exp(1));
ddiag = zeros(N,1);
ddiag(1) = 0;
ddiag(end) = alpha - h/2 * x(end);
for ii = 2 : N-1
    ddiag(ii) = h^2 * x(ii);
end
%===================================
% Matrix A
%===================================
A = zeros(N,N);
for i0 = 1 : N
    for j0 = 1 : N
        if(i0==j0+1)
            A(i0,j0) = adiag(i0);
        end
        if(j0==i0+1)
            A(i0,j0) = cdiag(i0);
        end
        if(i0==j0)
            A(i0,j0) = bdiag(i0);
        end
    end
end
%===================================
% Progonka method
%===================================
y = Progon(A,ddiag)';
%===================================
% Plot
%===================================
figure(1)
plot(x,y,'b','LineWidth',3)
hold on
grid on
plot(x,ur,'g--','LineWidth',3)
xlabel('\bf{x}')
ylabel('\bf{Solution}')
legend('\bf{Approximated y}', ...
    '\bf{Exact u}')
%===================================
% Error
%===================================
err = abs(y-ur);
display(['Maximal error: ',num2str(max(err))])
time = toc;
display(['Elapsed time: ',num2str(time)])