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

figure()
scatter3(term1, term2, y, 'filled')
hold on
x1fit = min(term1):.01:max(term1);
x2fit = min(term2):.01:max(term2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
hold off

dt_critical = [(b(1)/(2*b(2))*dx(1)^2)^(1/3);
               (b(1)/(2*b(2))*dx(2)^2)^(1/3);
               (b(1)/(2*b(2))*dx(3)^2)^(1/3);
               (b(1)/(2*b(2))*dx(4)^2)^(1/3)]

           
%% PART II: check the convergence rate for dt:

dx = M(2:5,1);
dt = dx.*dx;
err = M(2:5,6);
figure()
error_loglog(dt, err)
error_table(dt, err)
 

%% PART III: check the convergence rate for dx:
M2 = csvread('check_dx.csv');
dx = M2(:,1);
err1 = M2(:,2); % l1 norm
err2 = M2(:,3); % l2 norm
err3 = M2(:,4); % inf norm
figure()
error_loglog(dx, err3)
error_table(dx, err3)



