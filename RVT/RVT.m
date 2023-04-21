function [U,V,X,Y,sigx,sigy,sigxy,alfa,sigz,tmax] = ...
    RVT2(a,b,pr,rg,ts,uinc,umin,umax,vinc,vmin,vmax)
%RVT2 calculates the stress due to distant tectonic forces and gravity in
%   symmetric ridges an valleys
% [U,V,X,Y,sigx,sigy,sigxy,alfa,sigz,tmax] = ...
%    RVT2(a,b,pr,rg,ts,uinc,umin,umax,vinc,vmin,vmax)
%
% ************************************************************************
%                                  RVT.FOR
%
% This program was developed and written by W.Z. Savage and P.S. Powers for
% the elastic solution of tectonic and gravity stresses in isolated
% symmetric ridges and valleys.
%
% ************************************************************************
%
% RVT.FOR modified by Nate to run efficiently and accurately in MATLAB
%
% Inputs:
% -   a, distance from the center lone fo ridge or valley to inflection
%     point on the surface
% -   b, maximum ridge height (positive value) or valley depth (negative
%     value)
% -   pr, Poisson's ratio
% -   rg, density times acceleration due to gravity
% -   ts, far-field tectonic stress; compression (negative value), tension
%     (positive value)
% -   uinc, increment value of the u-coordinate
% -   umin, minimum u-coordinate
% -   umax, maximum u-coordinate
% -   vinc, increment value of the v-coordinate
% -   vmin, minimum v-coordinate
% -   vmax, maximum v-coordinate
%
% Outputs:
% -   U, horizontal distance
% -   V, vertical distance
% -   X, conformally mapped horizontal distance
% -   Y, conformally mapped vertical distance
% -   sigx, xx component of stress @ (X,Y)
% -   sigy, yy component of stress @ (X,Y)
% -   sigxy, xy component of stress @ (X,Y)
% -   alfa, orientation of principal stress @ (X,Y)
% -   sigz, zz component of stress @ (X,Y) (i.e. out of screen)
% -   tmax, maximum shear stress @ (X,Y)
%
% Last Edit: 2.7.23, Nate

% rad is the small radius from the point u = -ia
rad = 1e-4;

so = rg*b;
po = complex(0,1);
ai = po*a;

% Phi @ w = -ia for gravity solution
phia = -so*((4*a+b)*(1-pr)+b)/(8*(1-pr)*(2*a+b));

% Phi @ w = -ia for tectonic solution
phiat = -(ts*b)/(4*(2*a+b));

% d2phia is the 2nd derivative of phi at w = -ia for tectonic sol.
d2phia = -(ts*b*(4*a+b)*(b-12*a));
d2phia = d2phia/(2*a*(2*a+b)*(4*a+b)^3);

ulist = umin:uinc:umax;
vlist = vmin:vinc:vmax;

[U, V] = meshgrid(ulist, vlist);

w = complex(U,V);
r = sqrt(U.^2+(V+a).^2);

% conformal map
z = w+(a*b)./(w-ai);

X = real(z);
Y = imag(z);
dz = ((w-ai).^2-a*b)./((w-ai).^2);

% aw is a(w)
aw1 = po*(4*a+b)./(8*(w-ai));
aw2 = a*b*(w-3*ai)./(8*(1-pr)*(w-ai).^3);
aw = -aw1-aw2;

% phi is phi(w) for the gravity solution
phi = -aw*so./dz+a*b*phia./(dz.*(w-ai).^2);

% awt is A0(w)
awt = -(a*b*ts)./(2*(w-ai).^2);

% phit is phi0(w)
phit = -awt./dz+a*b*phiat./(dz.*(w-ai).^2);

% sum(sigmaxx + sigmayy) for tectonic (ssumt) and gravity (ssum)
ssumt = 4*real(phit)+ts;
ssum = ssumt+4*real(phi)+rg*Y/(1-pr);

% 1st derivative of phi(w)
d1aw1 = po*(4*a+b)./(8*(w-ai).^2);
d1aw2 = 2*a*b./(8*(1-pr)*(w-ai).^3);
d1aw3 = 6*po*a^2*b./(8*(1-pr)*(w-ai).^4);
d1aw = d1aw1+d1aw2-d1aw3;
d1phi = -so*d1aw./dz-(2*a*b*(phi+phia))./(dz.*((w-ai).^3));

% 2nd derivative of phi(w) to be used @ w = -ia
d2phi1 = 2*phi./((w-ai).^2);
d2phi2 = 4*d1phi./(w-ai);
d2phi3 = so*po*a^2*b./(2*(1-pr)*((w-ai).^5));
d2phi = -(d2phi1+d2phi2+d2phi3)./dz;

% B(w)
bw1 = po*(4*a+b)./(8*(w-ai));
bw2 = (1-2*pr)*a*b*(w-3*ai)./(8*(1-pr)*(w-ai).^3);
bw = -so*(bw1+bw2);

psi1 = w.*d1phi+bw+phi;
d1phi1 = -(a*b*ts*(4*a+b)*(w-ai));
d1phi2 = 2*(2*a+b)*(((w-ai).^2-a*b).^2);
d1phit = d1phi1./d1phi2;
bwt = -awt;

% part of (sigmayy - sigmaxx + 2isigmaxy)
psi1t = w.*d1phit+bwt+phit;

% Test on closeness of w to -ia.
% If w is near -ia, the Taylor expansion about -ia is used
psi2 = zeros(size(U));
psi2t = zeros(size(U));

% Taylor's expansion for gravity and tectonic stress at w = -ia
psi2(r < rad) = 0.5*a*b*d2phi(r < rad);
psi2t(r < rad) = 0.5*a*b*d2phia;
% Else
psi2(~(r < rad)) = a*b*d1phi(~(r < rad))./(w(~(r < rad))+ai)-...
    a*b*(phi(~(r < rad))-phia)./((w(~(r < rad))+ai).^2);
psi2t(~(r < rad)) = a*b*d1phit(~(r < rad))./(w(~(r < rad))+ai)-...
    a*b*(phit(~(r < rad))-phiat)./((w(~(r < rad))+ai).^2);

psi = -(psi1+psi2)./dz;
psit = -(psi1t+psi2t)./dz;
zbar = conj(z);

str1 = 2*(zbar.*d1phi./dz+psi);
str1t = 2*(zbar.*d1phit./dz+psit);

% Equations in this block form the sums and differences
dift = real(str1t) - ts;
dif = dift + real(str1) + rg*Y*(1-2*pr)/(1-pr);
sigxy = imag(str1) + imag(str1t);
sigxy = sigxy/2;
sigx = (ssum - dif)/2;
sigy = (ssum + dif)/2;

% Principal stress directions
alfa = 0.5*atan2((2*sigxy),(sigx-sigy));
alfa = (180/pi)*alfa;

% Principal stress magnitudes
sig1 = ssum/2+sqrt((dif/2).^2+sigxy.^2);
sig2 = ssum/2-sqrt((dif/2).^2+sigxy.^2);

% Maximum shear stress
tmax = (sig1 - sig2)/2;

% Out of plane stress
sigz = pr*(sigx+sigy);

end