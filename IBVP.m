tic
%===================================
% Boundary value problem (BVP) for
% ODE of 2 order
% u'' = u' + 2*y + cos(x)
% (1) u(0) = -0.3
% (2) u(pi/2) = -0.1
% u_{exact} = -1/10 * (sin(x) + 3*cos(x))
%===================================
clear all
clc
%===================================
% Grid
%===================================
a = 0; b = pi/2; h = 0.01;
x = a : h : b;
N = length(x);
y = zeros(1,N);
%===================================
% Exact solution
%===================================
ur = zeros(1,N);
for i = 1 : N
    ur(i) = -1/10 * (sin(x(i)) + ...
        3*cos(x(i)));
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
% Diagonals
%===================================
adiag = zeros(1,N);
bdiag = zeros(1,N);
cdiag = zeros(1,N);
ddiag = zeros(N,1);
%===================================
% adiag
%===================================
adiag(1) = 0;
adiag(end) = 0;
for ii = 2 : N-1
    adiag(ii) = 1 + h/2;
end
%===================================
% bdiag
%===================================
bdiag(1) = 1;
bdiag(end) = 1;
for ii = 2 : N-1
    bdiag(ii) = -2*(1+h^2);
end
%===================================
% cdiag
%===================================
cdiag(1) = 0;
cdiag(end) = 0;
for ii = 2 : N-1
    cdiag(ii) = 1 - h/2;
end
%===================================
% ddiag
%===================================
ddiag(1) = ur(1);
ddiag(end) = ur(end);
for ii = 2 : N-1
    ddiag(ii) = h^2 * cos(x(ii));
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
legend('\bf{Approximated}','\bf{Exact}')
%===================================
% Error
%===================================
err = abs(y-ur);
display(['Maximal error: ',num2str(max(err))])
time = toc;
display(['Elapsed time: ',num2str(time)])