
function y = Interpolate1D2D(xn,Vn,cl,xi)

Nxi = length(xi);
y = zeros(1,Nxi);

Ne = size(cl,1);

for i = 1:Nxi

    % find corresponding element
    for j = 1:Ne
        x1e = xn(cl(j,1));
        x2e = xn(cl(j,2));    
        x3e = xn(cl(j,3));
        if (xi(i) >= x1e) && (xi(i) <= x2e)
            break;
        end
    end

    % calculate the coordinate of xi in reference element
    zeta = 2*(xi(i) - x3e)/(x2e - x1e);

    % calculate the value of V
    N1 = zeta*(zeta - 1)/2;
    N2 = zeta*(zeta + 1)/2;
    N3 = (zeta+1)*(1-zeta);
    
    y(i) = Vn(cl(j,1))*N1 + Vn(cl(j,2))*N2 + Vn(cl(j,3))*N3;
end

end