clc;
clear;
close all;

syms zeta e le rho0

N1 = (1-zeta)/2;
N2 = (1+zeta)/2;

% K, stiffness matrix
(2/le)*int(diff(N1,zeta)*e*diff(N1,zeta),zeta,-1,1)
(2/le)*int(diff(N1,zeta)*e*diff(N2,zeta),zeta,-1,1)
(2/le)*int(diff(N2,zeta)*e*diff(N1,zeta),zeta,-1,1)
(2/le)*int(diff(N2,zeta)*e*diff(N2,zeta),zeta,-1,1)

% f, load vector
(-le*rho0/2)*int(N1,zeta,-1,1)
(-le*rho0/2)*int(N2,zeta,-1,1)
