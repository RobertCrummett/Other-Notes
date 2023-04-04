% Plot the analytic solution to Poissons equation
% Capacitor plates normal to x direction, d apart

clc; clear; close all;

% inputs
e0 = 8.85e-12; % F m^-1, electric perm of free space
er = 1; % relative electric perm
e = e0*er;
V0 = 1; % volts, voltage of plate at x=0
d = 8; % cm, distance between plates
rho0 = 1e-8; % C m^-3

% body
d = d/100; % converting units of d from cm to m
x = linspace(0,d,1000);
Vx = (d^2*rho0/(12*er*e0))*(1-x/d).^4 + (rho0*d/(12*er*e0) - V0/d)*x + (V0 - d^2*rho0/(12*er*e0));

% plot
figure()
plot(x,Vx,'b--')
xlabel("Distance (m)")
ylabel("Electric Potential (V)")