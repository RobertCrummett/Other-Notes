
function y = Interpolate1D1D(xn,Vn,cl,xi)

Nxi = length(xi);
y = zeros(1,Nxi);

Ne = size(cl,1);

for i = 1:Nxi

    % find corresponding element
    for j = 1:Ne
        x1e = xn(cl(j,1));
        x2e = xn(cl(j,2));    
        if (xi(i) >= x1e) && (xi(i) <= x2e)
            break;
        end
    end

    % calculate the coordinate of xi in reference element
    zeta = 2*(xi(i) - x1e)/(x2e - x1e) - 1;

    % calculate the value of V
    N1 = (1-zeta)/2;
    N2 = (1+zeta)/2;
    
    y(i) = Vn(cl(j,1))*N1 + Vn(cl(j,2))*N2;
end

end