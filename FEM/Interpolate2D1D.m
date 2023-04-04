function vi = Interpolate2D1D(p,cl,Vn,dp)

TR = triangulation(cl',p');
tid = pointLocation(TR,dp);

Np = length(tid);
vi = zeros(1,Np);


for i = 1:Np
    
    % calculate the coordinate of dp in reference coordinate system
    x = dp(i,1);
    y = dp(i,2);
    x1 = p(1,cl(1,tid(i)));
    y1 = p(2,cl(1,tid(i)));
    x21 = p(1,cl(2,tid(i))) - x1;
    x31 = p(1,cl(3,tid(i))) - x1;
    y21 = p(2,cl(2,tid(i))) - y1;
    y31 = p(2,cl(3,tid(i))) - y1;

    det = x21*y31 - x31*y21;
    zeta = (y31*(x-x1) - x31*(y-y1))/det;
    eta = (-y21*(x-x1) + x21*(y-y1))/det;

    % calculate the value of V
    N1 = 1-zeta-eta;
    N2 = zeta;
    N3 = eta;
    
    vi(i) = Vn(cl(1,tid(i)))*N1 + Vn(cl(2,tid(i)))*N2 + Vn(cl(3,tid(i)))*N3;
end

end