clc;
clear;
close all;

V0 = 1;
W = 1;
H = 1;

Nx = 10;
Ny = 10;

dx = W/(Nx-1);
dy = H/(Ny-1);

% allocate memory for mesh nodes & connectivity list
p = zeros(2,Nx*Ny);

index = 0;
for i = 1:Ny
    y = ((i-1)*H)/(Ny-1);
    for j = 1:Nx
        x = ((j-1)*W)/(Nx-1);
        index = index + 1;
        p([1,2],index) = [x,y];
    end
end

cl = delaunay(p')'; % good for convex domain

patch('faces',cl','vertices',p','facecolor','c','edgecolor','k')
axis off equal

% Analytic solution
Vex = zeros(Nx*Ny,1);

for k = 1:100
    Vex = Vex + (sin((2*k-1)*pi*p(1,:)'/W).*sinh((2*k-1)*pi*p(2,:)'/W))...
        ./((2*k-1)*sinh((2*k-1)*pi*H/W));
    % sinh -> inf as k -> 110, stopped infinite summation at k = 100
end
Vex = Vex*4*V0/pi;

figure;
patch('faces',cl','vertices',p','edgecolor','none','facecolor','interp','facevertexcdata',Vex)
colorbar
colormap('jet')
title("Analytic Solution")


% allocate memory for stiffness matrix and load vector
% using sparse matrix
Ne = size(cl,2);
Nn = Nx*Ny;
IndexI = zeros(9,Ne);
IndexJ = zeros(9,Ne);
valK = zeros(9,Ne);
f = zeros(Nn,1);
Ae = dx*dy/2; % area of every element is the same

x13 = p(1,cl(1,:)) - p(1,cl(3,:));
x31 = -x13;
x23 = p(1,cl(2,:)) - p(1,cl(3,:));
x32 = -x23;
x12 = p(1,cl(1,:)) - p(1,cl(2,:));
x21 = -x12;
y12 = p(2,cl(1,:)) - p(2,cl(2,:));
y21 = -y12;
y23 = p(2,cl(2,:)) - p(2,cl(3,:));
y32 = -y23;
y31 = p(2,cl(3,:)) - p(2,cl(1,:));
y13 = -y31;

% M11
IndexI(1,:) = cl(1,:);
IndexJ(1,:) = cl(1,:);
valK(1,:) = (-1./(4*Ae)).*(y23.^2 + x32.^2);

% M21
IndexI(2,:) = cl(2,:);
IndexJ(2,:) = cl(1,:);
valK(2,:) = (-1./(4*Ae)).*(y23.*y31 + x32.*x13);

% M12
IndexI(4,:) = cl(1,:);
IndexJ(4,:) = cl(2,:);
valK(4,:) = valK(2,:);

% M31
IndexI(3,:) = cl(3,:);
IndexJ(3,:) = cl(1,:);
valK(3,:) = (-1./(4*Ae)).*(y23.*y12 + x32.*x21);

% M13
IndexI(7,:) = cl(1,:);
IndexJ(7,:) = cl(3,:);
valK(7,:) = valK(3,:);

% M22
IndexI(5,:) = cl(2,:);
IndexJ(5,:) = cl(2,:);
valK(5,:) = (-1./(4*Ae)).*(y31.^2 + x13.^2);

% M33
IndexI(9,:) = cl(3,:);
IndexJ(9,:) = cl(3,:);
valK(9,:) = (-1./(4*Ae)).*(y12.^2 + x21.^2);

% M32
IndexI(6,:) = cl(3,:);
IndexJ(6,:) = cl(2,:);
valK(6,:) = (-1./(4*Ae)).*(y31.*y12 + x13.*x21);

% M23
IndexI(8,:) = cl(2,:);
IndexJ(8,:) = cl(3,:);
valK(8,:) = valK(6,:);

% stiffness matrix
K = sparse(IndexI,IndexJ,valK);

% imposing Dirichlet boundary conditions
% top points
index = Nx*(Ny-1)+1:Nx*Ny;
K(index,:) = 0;
for i = index
    K(i,i) = 1;
end
% f(index) = V0*sin(pi*p(1,index));
f(index) = V0;

% bottom index
index = abs(p(2,:) - 0) < 1e-6;
index = bitor(index,abs(p(1,:) - 0) < 1e-6);
index = bitor(index,abs(p(1,:) - W) < 1e-6);
tmp = 1:Nn;
index = tmp(index);
K(index,:) = 0;
for i = index
    K(i,i) = 1;
end
f(index) = 0; % should this be there?

% solving system of equations
V = K\f;

figure;
patch('faces',cl','vertices',p','edgecolor','none','facecolor','interp',...
    'facevertexcdata',V,'edgecolor','none')
colorbar
colormap('jet')
title("FEM Solution")

figure;
index = 11:Nx:Nx*Ny+1;
plot(linspace(0,1,length(index)),Vex(fliplr(index)),'r')
hold on
plot(linspace(0,1,length(index)),V(fliplr(index)),'b')
hold off
