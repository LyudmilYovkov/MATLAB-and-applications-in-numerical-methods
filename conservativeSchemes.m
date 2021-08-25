tic
%===================================
% Conservative difference scheme
% for ODE of 2 order with Dirichlet
% boundary conditions
% (k(x)u'(x))' - q(x)u(x) = -f(x)
% u_{exact} = sin(x)-x^2
% k(x) = x^2+1
% q(x) = x
% f(x) = -x^3 + 6*x^2 + (x^2+x+1)*sin(x) - 2*x*cos(x) + 2
% u(0) = u_{exact}(0) = 0
% u(1) = u_{exact}(1) = sin(1)-1 = -0.1585
%===================================
clear;
clc;
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
ur = zeros(1,N);
for i = 1 : N
    ur(i) = sin(x(i)) - x(i)^2;
end
%===================================
% Left boundary condition
%===================================
y(1) = ur(1);
%===================================
% Right boundary condition
%===================================
y(N) = ur(N);
%===================================
% Function k(x)
%===================================
kF = inline('x.^2+1','x');
%===================================
% Function q(x)
%===================================
qF = inline('x','x');
%===================================
% Function f(x)
%===================================
fF = inline('-x.^3+6*x.^2+(x.^2+x+1).*sin(x)-2*x.*cos(x)+2','x');
%===================================
% adiag
%===================================
adiag = zeros(1,N);
adiag(1) = 0;
adiag(end) = 0;
for ii = 2 : N-1
    adiag(ii) = kF(x(ii)-h/2);
end
%===================================
% bdiag
%===================================
bdiag = zeros(1,N);
bdiag(1) = 1;
bdiag(end) = 1;
for ii = 2 : N-1
    bdiag(ii) = -(kF(x(ii)+h/2) + ...
        kF(x(ii)-h/2) + h^2 * qF(x(ii)));
end
%===================================
% cdiag
%===================================
cdiag = zeros(1,N);
cdiag(1) = 0;
cdiag(end) = 0;
for ii = 2 : N-1
    cdiag(ii) = kF(x(ii)+h/2);
end
%===================================
% ddiag
%===================================
ddiag = zeros(N,1);
ddiag(1) = y(1);
ddiag(end) = y(end);
for ii = 2 : N-1
    ddiag(ii) = -h^2 * fF(x(ii));
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
% Solving by Progonka
%===================================
y = Progon(A,ddiag)';
%===================================
% Plot
%===================================
figure(1)
plot(x,y,'b','LineWidth',3)
hold on
grid on
plot(x,ur,'cyan--','LineWidth',3)
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