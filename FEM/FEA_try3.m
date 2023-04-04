clear; clf; close all
% Multi-domain problem
AnalyticSolution;
Ne = 100;
err = [];
for Nn = Ne+1 % number of mesh nodes

    xn = linspace(0,d,Nn); % x coordinate of each mesh node
    yn = zeros(1,Nn);
    p = [xn' yn'];

    % calculate mesh length
    Ne = Nn - 1; % number of elements
    cl = zeros(Ne,2); % connectivity list
    for i = 1:Ne
        for j = 1:2
            cl(i,j) = i + j - 1;
            % records corresponding nodes of each element
        end
    end

    % calculate e
    e = zeros(1,Ne);
    for i = 1:Ne
        xc = (xn(cl(i,2)) + xn(cl(i,1)))/2;
        if xc < d/2
            e(i) = 1*8.85e-12;
        else
            e(i) = 12*8.85e-12;
        end
    end

    % calculate mesh length
    le = zeros(1,Ne);
    for i = 1:Ne
        le(i) = xn(cl(i,2)) - xn(cl(i,1));
    end

    hold on
    patch('Faces',cl,'Vertices',p)
    % plot nodes
    plot(xn,yn,'yo',MarkerSize=15)

    % allocate memory for stiffness matrix and load vector
    K = zeros(Nn,Nn); % later we will use sparse
    f = zeros(Nn,1);

    % global assembly
    for i = 1:Ne
        % stiffness matrix
        K(cl(i,1),cl(i,1)) = K(cl(i,1),cl(i,1)) + (e(i)/le(i));
        K(cl(i,1),cl(i,2)) = K(cl(i,1),cl(i,2)) - (e(i)/le(i));
        K(cl(i,2),cl(i,1)) = K(cl(i,2),cl(i,1)) - (e(i)/le(i));
        K(cl(i,2),cl(i,2)) = K(cl(i,2),cl(i,2)) + (e(i)/le(i));
        % force vector
        f(cl(i,1)) = f(cl(i,1)) - le(i)*rho0/2;
        f(cl(i,2)) = f(cl(i,2)) - le(i)*rho0/2;
    end

    % imposing Dirichlet boundary conditions
    % left point
    K(1,1) = 1;
    for j = 2:Nn
        K(1,j) = 0;
    end
    f(1) = V0;
    % right point
    K(Nn,Nn) = 1;
    for j = 1:Nn-1
        K(Nn,j) = 0;
    end
    f(Nn) = 0;

    % solving system of equations
    V = K\f;

    Vfe = Interpolate1D1D(xn,V,cl,x);

    plot(x,Vfe,'color','red')

end