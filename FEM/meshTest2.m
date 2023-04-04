clc;
clear;
close all;

p = [0 2 2 1 1 0; 0 0 2 2 1 1]; % non convex region

cl = delaunay(p')';

xc = (p(1,cl(1,:))+p(1,cl(2,:))+p(1,cl(3,:)))/3;
yc = (p(2,cl(1,:))+p(2,cl(2,:))+p(2,cl(3,:)))/3;

in = inpolygon(xc,yc,p(1,:),p(2,:)); % center points within region
cl = cl(:,in);

hold on
patch('faces',cl','vertices',p','facecolor','c','edgecolor','k')
plot(p(1,:),p(2,:),'o','color','k','MarkerFaceColor','w','MarkerSize',7)
axis off equal

plot(xc,yc,'*')

