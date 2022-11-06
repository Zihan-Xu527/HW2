clc; clear all; close all;

%% PART I: regress err vs dx^2/dt and dt^2
M=csvread('inf_errors.csv');
dx = [M(2:5,1);M(2:5,1);M(2:5,1);M(2:5,1)];
dt = [M(2:5,1).* 2; 
      M(2:5,1); 
      M(2:5,1)./5; 
      M(2:5,1)./10];
y = [M(2:5,2);M(2:5,3);M(2:5,4);M(2:5,5)]; % errors  
dx2 = dx.*dx;


% regress errors in terms of dx^2/dt and dt^2
term1 = dx2./dt;
term2 = dt.*dt;
X = [ones(size(term1)) term1 term2];
b = regress(y, X)


dt_critical = [(b(2)/(2*b(3))*dx(1)^2)^(1/3);
               (b(2)/(2*b(3))*dx(2)^2)^(1/3);
               (b(2)/(2*b(3))*dx(3)^2)^(1/3);
               (b(2)/(2*b(3))*dx(4)^2)^(1/3)]


figure()
% regression estimation
t2=linspace(dt(2),dt(14));
y2=b(1)+b(2)*dx(2)^2*t2.^(-1)+t2.^2; 
subplot(1,3,1)
plot(t2,y2,'r--');
hold on
% critical point based on regression
y_critical =b(1)+b(2)*dx(2)^2*dt_critical(2)^(-1)+dt_critical(2)^2; 
scatter(dt_critical(2),y_critical,'r');
hold on
% true err (input data)
t2=[dt(2);dt(6);dt(10);dt(14)];
plot(t2,M(3,2:5),'bx-') 

xlabel('\Deltat')
ylabel('||e||_{\infty}')
legend({'estimated error','true error','critical \Deltat'})
title("Errors vs \Deltat with \Deltax=1/16")
hold on

% regression estimation
t3=linspace(dt(3),dt(15));
y3=b(1)+t3.^(-1)*b(2)*dx(3)^2+t3.^2;
subplot(1,3,2)
plot(t3,y3,'r--');
hold on
% critical point based on regression
y_critical =b(1)+b(2)*dx(3)^2*dt_critical(3)^(-1)+dt_critical(3)^2;
scatter(dt_critical(3),y_critical,'r')
hold on
% true err (input data)
t3=[dt(3);dt(7);dt(11);dt(15)];
plot(t3,M(4,2:5),'bx-')

xlabel('\Deltat')
ylabel('||e||_{\infty}')
legend({'estimated error','true error','critical \Deltat'})
title("Errors vs \Deltat with \Deltax=1/32")

% regression estimation
t4=linspace(dt(4),dt(16));
y4=b(1)+t4.^(-1)*b(2)*dx(4)^2+t4.^2;
subplot(1,3,3)
plot(t4,y4,'r--');
hold on
% critical point based on regression
y_critical =b(1)+b(2)*dx(4)^2*dt_critical(4)^(-1)+dt_critical(4)^2;
scatter(dt_critical(3),y_critical,'r')
hold on
% true err (input data)
t4=[dt(4);dt(8);dt(12);dt(16)];
plot(t4,M(5,2:5),'bx-')

xlabel('\Deltat')
ylabel('||e||_{\infty}')
legend({'estimated error','true error','critical \Deltat'})
title("Errors vs \Deltat with \Deltax=1/64")
hold off

           
% %% PART II: check the convergence rate for dt:
% 
% dx = M(2:5,1);
% dt = dx.*dx;
% err = M(2:5,6);
% figure()
% error_loglog(dt, err)
% error_table(dt, err)
%  
% 
% %% PART III: check the convergence rate for dx:
% M2 = csvread('check_dx.csv');
% dx = M2(:,1);
% err1 = M2(:,2); % l1 norm
% err2 = M2(:,3); % l2 norm
% err3 = M2(:,4); % inf norm
% figure()
% error_loglog(dx, err3)
% error_table(dx, err3)



