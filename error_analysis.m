clc; clear all;
M=csvread('inf_errors.csv');
M
dx=M(2:5,1);
dt = dx/5;
y=M(2:5,4).*dt;
% loglog(x,y);
% polyfit(log(x),log(y),1)
error_loglog(dx,y)
error_table(dx,y)