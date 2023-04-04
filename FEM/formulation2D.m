clc;
clear;
close all;

syms alpha_x alpha_y e z x1 y1 x21 y21 x31 y31 x23 y23

x = x1 + x21*z + x31*e;
y = y1 + y21*z + y31*e;

% Jacobian matrix
J = [diff(x,z) diff(y,z); diff(x,e) diff(y,e)];
invJ = inv(J);

N = [1-z-e z e];

rNrX = invJ(1,:)*[diff(N,z);diff(N,e)];
rNrY = invJ(2,:)*[diff(N,z);diff(N,e)];

for j = 1:3
    for i = 1:3
        tmp = -int(int((alpha_x*rNrX(i)*rNrX(j) + ...
            alpha_y*rNrY(i)*rNrY(j))*det(J),z,0,1-e),e,0,1);
        tmp = subs(tmp,x21,x23+x31);
        tmp = subs(tmp,y21,y23+y31);
        pretty(simplify(tmp))
    end
end
