clc;
clear;
close all;

x = [1 3 5];
y = [1 1 4];

hold on
subplot(2,2,2)
plot(x,y,'ok')
patch('faces',[1 2 3],'Vertices',[x;y]','facecolor','c','edgecolor','k')
axis square

x21 = x(2) - x(1);
x31 = x(3) - x(1);
y21 = y(2) - y(1);
y31 = y(3) - y(1);

e = [0 1 0];
z = [0 0 1];

subplot(2,2,1)
plot(e,z,'ok')
patch('faces',[1 2 3],'Vertices',[e;z]','facecolor','m','edgecolor','k')
axis square

xm = x(1) + x21*z + x31*e;
ym = y(1) + y21*z + y31*e;

subplot(2,2,4)
plot(xm,ym,'ok')
patch('faces',[1 2 3],'Vertices',[x;y]','facecolor','r','edgecolor','k')
axis square
