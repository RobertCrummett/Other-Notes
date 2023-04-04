function y = CalculateE1D2D(xn,Vn,cl,xi)
Nxi = length(xi);
y = zeros(1,Nxi);

Ne = size(cl,1);

for i = 1:Nxi
    % find the corresponding element
    for j = 1:Ne
        x1e = xn(cl(j,1));
        x2e = xn(cl(j,2));
        x3e = xn(cl(j,3));

        if (xi(i) >= x1e) && (xi(i) <= x2e)
            break
        end
    end

    V1e = Vn(cl(j,1));
    V2e = Vn(cl(j,2));
    V3e = Vn(cl(j,3));

    zeta = 2*(xi(i) - x3e)/(x2e - x1e);

    y(i) = (-2/(x2e - x1e))*(V1e*(zeta-1/2) + V2e*(zeta+1/2) ...
        + V3e*(-2*zeta));
end
end