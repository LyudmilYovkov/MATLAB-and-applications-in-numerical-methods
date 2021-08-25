%===================================
% Nonlinear one-dimensional scheme
% for the heat equation 
% u^2 * u_{t} = (u * u_{x})_{x} + u_{zz} + f
% 0 <= x <= 1, 0 <= z <= 1, 0 <= t <= 1
% u(t=0) = sin(x) + cos(z)
% u(x=0) = cos(z) + t^4
% u(x=1) = sin(1) + cos(z) + t^4
% u(z=0) = sin(x) + 1 + t^4
% u(z=1) = sin(x) + cos(1) + t^4
% f = 4 * t^3 * u^2 - cos(x)^2 + u * sin(x) + cos(z)
% u_{exact} = sin(x) + cos(z) + t^4
%===================================
clear;
clc;
close;
tic;
%===================================
% Grid
%===================================
ax = 0; bx = 1; hx = (bx-ax)/60;
x = ax : hx : bx;
N = length(x);
az = 0; bz = 1; hz = (bz-az)/50;
z = az : hz : bz;
M = length(z);
t0 = 0; tend = 1; tau = (tend-t0)/600;
t = t0 : tau : tend;
T = length(t);
%===================================
% Exact solution
%===================================
ur = zeros(N,M,T);
for i = 1 : N
    for j = 1 : M
        for n = 1 : T
            ur(i,j,n) = sin(x(i)) + ...
                cos(z(j)) + t(n)^4;
        end
    end
end
%===================================
% Approximating solution 
%===================================
u = zeros(N,M,T);
%===================================
% Initial condition
%===================================
for i = 1 : N
    for j = 1 : M
        u(i,j,1) = sin(x(i)) + cos(z(j));
    end
end
%===================================
% Locally 1D scheme
%===================================
for n = 1 : T-1
    sol = zeros(N,M);
    %===================================
    % 1D problem for x
    %===================================
    for j = 1 : M
        %===================================
        % adiag
        %===================================
        adiag = zeros(1,N);
        adiag(1) = 0;
        adiag(end) = 0;
        for ii = 2 : N-1
            adiag(ii) = -tau * 0.5 * (u(ii-1,j,n) + u(ii,j,n)) / ...
                ((u(ii,j,n)^2) * (hx^2));
        end
        %===================================
        % bdiag
        %===================================
        bdiag = zeros(1,N);
        bdiag(1) = 1;
        bdiag(end) = 1;
        for ii = 2 : N-1
            bdiag(ii) = 1 + tau/((u(ii,j,n)^2) * (hx^2)) * ...
                0.5 * (u(ii-1,j,n) + 2 * u(ii,j,n) + u(ii+1,j,n));
        end
        %===================================
        % cdiag
        %===================================
        cdiag = zeros(1,N);
        cdiag(1) = 0;
        cdiag(end) = 0;
        for ii = 2 : N-1
            cdiag(ii) = -tau * 0.5 * (u(ii+1,j,n) + u(ii,j,n)) / ...
                ((u(ii,j,n)^2) * (hx^2));
        end
        %===================================
        % ddiag
        %===================================
        ddiag = zeros(N,1);
        ddiag(1) = cos(z(j)) + (t(n)+tau/2)^4;
        ddiag(end) = sin(1) + cos(z(j)) + (t(n)+tau/2)^4;
        for ii = 2 : N-1
            ddiag(ii) = u(ii,j,n) + tau/(2*(u(ii,j,n)^2)) * ...
                (4 * (t(n)+tau/2)^3 * u(ii,j,n)^2 - cos(x(ii))^2 + ...
                u(ii,j,n) * sin(x(ii)) + cos(z(j)));
        end
        %===================================
        % Matrix A
        %===================================
        A = zeros(N,N);
        for i0 = 1 : N
            for j0 = 1 : N
                if(i0==j0+1)
                    A(i0,j0) = adiag(i0);
                elseif(j0==i0+1)
                    A(i0,j0) = cdiag(i0);
                elseif(i0==j0)
                    A(i0,j0) = bdiag(i0);
                end
            end
        end
        %===================================
        % Progonka for the middle solution u(n+1/2)
        %===================================
        sol(:,j) = Progon(A,ddiag);
    end
    %===================================
    % 1D problem for z
    %===================================
    for i = 1 : N
        %===================================
        % adiag
        %===================================
        adiag = zeros(1,M);
        adiag(1) = 0;
        adiag(end) = 0;
        for jj = 2 : M-1
            adiag(jj) = -tau/(((sol(i,jj))^2) * (hz^2));
        end
        %===================================
        % bdiag
        %===================================
        bdiag = zeros(1,M);
        bdiag(1) = 1;
        bdiag(end) = 1;
        for jj = 2 : M-1
            bdiag(jj) = 1 + 2*tau/(((sol(i,jj))^2) * (hz^2));
        end
        %===================================
        % cdiag
        %===================================
        cdiag = zeros(1,M);
        cdiag(1) = 0;
        cdiag(end) = 0;
        for jj = 2 : M-1
            cdiag(jj) = -tau/(((sol(i,jj))^2) * (hz^2));
        end
        %===================================
        % ddiag
        %===================================
        ddiag = zeros(M,1);
        ddiag(1) = sin(x(i)) + 1 + (t(n+1))^4;
        ddiag(end) = sin(x(i)) + cos(1) + (t(n+1))^4;
        for jj = 2 : M-1
            ddiag(jj) = sol(i,jj) + tau/(2 * (sol(i,jj))^2) * ...
                (4 * (t(n)+tau/2)^3 * (sol(i,jj))^2 - cos(x(i))^2 + ...
                sol(i,jj) * sin(x(i)) + cos(z(jj)));
        end
        %===================================
        % Matrix A
        %===================================
        A = zeros(M,M);
        for i0 = 1 : M
            for j0 = 1 : M
                if(i0==j0+1)
                    A(i0,j0) = adiag(i0);
                elseif(j0==i0+1)
                    A(i0,j0) = cdiag(i0);
                elseif(i0==j0)
                    A(i0,j0) = bdiag(i0);
                end
            end
        end
        %===================================
        % Progonka for u(n+1)
        %===================================
        u(i,:,n+1) = Progon(A,ddiag);
    end
end
%===================================
% Plot
%===================================
plot(x,ur(:,end-3,end-10),'b','LineWidth',3)
hold on, grid on
plot(x,u(:,end-3,end-10),'g--','LineWidth',3)
legend('\bf{Exact}','\bf{Approximating}')
%===================================
% Comparing with the exact solution
%===================================
err = abs(ur-u);
Max = max(max(max(err)));
relErr = max(max(max(err)))/max(max(max(u)))*100;
display(['tau = ',num2str(tau)])
display(['Maximal absolute error: ',num2str(Max)])
display(['Relative error: ',num2str(relErr)])
time = toc;
display(['Elapsed time: ',num2str(time)])
display('==================================')