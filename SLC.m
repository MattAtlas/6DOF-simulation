clear all;
close all;
clc;

h = 0.0001;          % time step
t = 0:h:2;          % time range
N = length(t);      % number of steps
Iy = 0.007875;

x_inner = zeros(3,N);       % theta theta_dot motor_response
x_inner(:,1) = [0;0;0];
x_outer = zeros(2,N);       % x x_dot
x_outer(:,1) = [-1;0];

u_inner = zeros(1,N);
u_outer = zeros(1,N);
max_F = [1.704 19.472 1.014]; % x z theta

A_inner = [0 1 0;0 0 max_F(3)/Iy;0 0 -20]; 
B_inner = [0 0 20]';
C_inner = [1 0 0;0 1 0];
D_inner = [0 0]';

sys1 = ss(A_inner,B_inner,C_inner,D_inner,0);
[K_inner] = lqr(sys1,diag([10 1 1]),0.1)

A_outer = [0 1;0 0];
B_outer = [0 9.81]';
C_outer = [1 0];
D_outer = 0;

sys2 = ss(A_outer,B_outer,C_outer,D_outer,0);
[K_outer] = lqr(sys2,diag([40 2]),0.01)
%K_outer = [10 2]

for n=1:N-1

    u_outer(n) = -K_outer*x_outer(:,n);
    
        if u_outer(n) >= 0.5;
            u_outer(n) = 0.5;
        end
        if u_outer(n) <= -0.5;
            u_outer(n) = -0.5;
        end

    u_inner(n) = -K_inner*(x_inner(:,n)-[u_outer(n);0;0]);

        if u_inner(n) >= max_F(3)
            u_inner(n) = max_F(3);
        end
        if u_inner(n) <= -max_F(3);
            u_inner(n) = -max_F(3);
        end

    x_inner(:,n+1) = x_inner(:,n)+(A_inner*x_inner(:,n) + B_inner*u_inner(n))*h;
    
    x_outer(:,n+1) = x_outer(:,n)+(A_outer*x_outer(:,n) + B_outer*u_outer(n))*h;
    
    u_outer(n+1) = x_inner(1,n+1);

    
    
end
figure(1);
hold on;
%[AX,H1,H2] = plotyy(t(2:N-1),x(1,2:N-1),t(2:N-1),u(2:N-1),'plot');
plot(t(1:N-1),x_inner(1,1:N-1),'b');
plot(t(1:N-1),x_inner(2,1:N-1),'r');
plot(t(1:N-1),x_inner(3,1:N-1),'g');
plot(t(1:N-1),u_inner(1:N-1),'k');
legend('theta (rad)','theta dot (rad/s)','motor (Nm)','input u (Nm)');
title('Inner Loop Step Characteristics');
xlabel('Time in s');

hold off;

figure(2);
hold on;
plot(t(1:N-1),x_outer(1,1:N-1),'b');
plot(t(1:N-1),x_outer(2,1:N-1),'r');
plot(t(1:N-1),u_outer(1:N-1),'k');
legend('x (m)','x dot (m/s)','theta ref (rad)');
xlabel('Time in s');
title('Outer Loop Step Characteristics');
% xlabel('Time in s');
% 
% set(get(AX(1),'Ylabel'),'String','Pos in x')
% set(get(AX(2),'Ylabel'),'String','Normalized input u')

eig(A_inner)
eig(A_inner-B_inner*K_inner)
