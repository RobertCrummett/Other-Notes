AnalyticSolution;
err = [];
Ne = 1; % number of elements
for Nn = 2*Ne+1 % number of mesh nodes

    xn = linspace(0,d,Nn); % x coordinate of each mesh node
    yn = zeros(1,Nn);
    p = [xn' yn'];

    % calculate mesh length
    cl = zeros(Ne,3); % connectivity list
    for i = 1:Ne
        cl(i,1) = 2*i-1;
        cl(i,2) = 2*i+1;
        cl(i,3) = 2*i;
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

        % Ke
        Ke = e/(3*le(i))*[7 1 -8; 1 7 -8; -8 -8 16];
        fe = (-le(i)*rho0/6)*[1;1;4];
        % stifness matix
        for j = 1:3
            for k = 1:3
                K(cl(i,j),cl(i,k)) = K(cl(i,j),cl(i,k)) + Ke(j,k);
            end
            f(cl(i,j)) = f(cl(i,j)) + fe(j);
        end
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

    % plot FEA
%     plot(xn,V,color='red')

    Vfe = Interpolate1D2D(xn,V,cl,x);

    plot(x,Vfe,'color','r')

    AreaError = sum(abs(Vx-Vfe))*100/abs(sum(Vx));
    err(end+1) = AreaError;
    
    L2Error = sqrt(sum((Vx-Vfe).^2*x(2)));

end