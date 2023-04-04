AnalyticSolution;
err = [];
for Nn = 300 % number of mesh nodes

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
    % using sparse matrix
    IndexI = zeros(4,Ne);
    IndexJ = zeros(4,Ne);
    valK = zeros(4,Ne);
    
    f = zeros(Nn,1);

    % global assembly
    for i = 1:Ne

        IndexI(:,i) = [cl(i,1);cl(i,1);cl(i,2);cl(i,2)];
        IndexJ(:,i) = [cl(i,1);cl(i,2);cl(i,1);cl(i,2)];
        valK(:,i) = [(e/le(i));-(e/le(i));-(e/le(i));(e/le(i))];

        % force vector
        f(cl(i,1)) = f(cl(i,1)) - le(i)*rho0/2;
        f(cl(i,2)) = f(cl(i,2)) - le(i)*rho0/2;
    end
    
    K = sparse(IndexI(:),IndexJ(:),valK(:));
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
    plot(xn,V,color = 'red')

    Vfe = Interpolate1D1D(xn,V,cl,x);

    AreaError = sum(abs(Vx - Vfe))*100/abs(sum(Vx));
    err(end+1) = AreaError;

    L2Error = sqrt(sum((Vx-Vfe).^2*x(2)));

end