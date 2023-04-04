function y = CalculateE1D1D(xn,Vn,cl,xi)
Nxi = length(xi);
y = zeros(1,Nxi);

Ne = size(cl,1);

for i = 1:Nxi
    % find the corresponding element
    for j = 1:Ne
        x1e = xn(cl(j,1));
        x2e = xn(cl(j,2));

        if (xi(i) >= x1e) && (xi(i) <= x2e)
            break
        end
    end

    zeta = 2*(xi(i) - x1e)/(x2e - x1e) - 1;

    y(i) = (Vn(cl(j,1)) - Vn(cl(j,2))) / (x2e - x1e);
end
end