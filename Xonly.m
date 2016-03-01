clear all;
close all;
clc;

h = 0.0001;          % time step
t = 0:h:2.5;          % time range
N = length(t);      % number of steps
Iy = 0.007875;
m = 1.01;
x_inner = zeros(3,N);       % theta theta_dot motor_response
x_inner(:,1) = [-1;0;0];
x_outer = zeros(2,N);       % x x_dot
x_outer(:,1) = [-1;0];

u_inner = zeros(1,N);
u_outer = zeros(1,N);
max_F = [1.704 19.472 1.014]; % x z theta

A_inner = [0 1 0;0 0 max_F(1)/m;0 0 -20]; 
B_inner = [0 0 20]';
C_inner = [1 0 0;0 1 0];
D_inner = [0 0]';

sys1 = ss(A_inner,B_inner,C_inner,D_inner,0);
[K_inner] = lqr(sys1,diag([300 10 1]),.1)

for n=1:N-1

    u_inner(n) = -K_inner*(x_inner(:,n));

        if u_inner(n) >= max_F(1)
            u_inner(n) = max_F(1);
        end
        if u_inner(n) <= -max_F(1);
            u_inner(n) = -max_F(1);
        end

    x_inner(:,n+1) = x_inner(:,n)+(A_inner*x_inner(:,n) + B_inner*u_inner(n))*h;
   
end
figure(1);
hold on;

plot(t(1:N-1),x_inner(1,1:N-1),'b');
plot(t(1:N-1),x_inner(2,1:N-1),'r');
plot(t(1:N-1),x_inner(3,1:N-1),'g');
plot(t(1:N-1),u_inner(1:N-1),'k');
legend('x (m)','x dot (m/s)','motor (N)','input u (N)');
title('6DOF Step Characteristics in x');
xlabel('Time (s)');

hold off;

print('-dpng','-r300','Xonly.png');
