function plotSphere(phi,the,field,lim)
% XYZ coordinates
[X,Y,Z] = sph2cart(phi,pi/2-the,1.002);

% Draw contours and sphere
M = contour(phi,pi/2-the,field); clf;
[XS,YS,ZS] = sphere(100);

% Draw sphere
surf(XS,YS,ZS,'FaceColor','white','EdgeColor','none')
hold on

% Draw contours
start = 1;
while start < size(M,2)
    count = M(2,start);
    [x,y,z] = sph2cart(M(1,(start+1):(start+count)),...
        M(2,(start+1):(start+count)),1.004);
    plot3(x,y,z,'k-','LineWidth',1.1)
    start = start + count + 1;
end

% Draw data
s = surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.75);

% Define RGB matrix to color data
CT = jet(256);
TN = round(rescale(double((field-min(1))/(lim(2) - lim(1))),1,size(CT,1)));
outpic = ind2rgb(TN,CT);
s.CData = outpic;

hold off
end