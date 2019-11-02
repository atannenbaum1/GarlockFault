function [F] = dipdikemodel(beta0,X)
% Magnetic field from a vertical dike using Telford Eq.3.44c

%X = dis

% I = inclination (degrees)
% Beta = strike (degrees) of dike anticlockwise from North magnetic
% x = Coordinates
% x0 = Dike location
% k = Susceptibility factor
% Fe = Earth's field
% w = Dike width
% z = Depth to top of dike
% L = vertical "length" of the dike

%%% Implement Eq 3.44c from Telford
%defined variables
%zetag = beta0(7);
I = 60*pi/180;
Fe = 48700; 
zeta = 90*pi/180;
declination = 12.4*pi/180;


% x0=a(1); %%% location
% k=a(2); %%% susceptibility
% z1=a(3); %%% depth to top
% off=a(4); %%% offset
% w=a(5); %%% dike width

x0=beta0(1);
k=beta0(2);
z=beta0(3);
w=beta0(4);
off=beta0(5);
strikeg = beta0(6);
strike = strikeg*pi/180; %measured from map
L=25000;
beta0;
%%% Bring angles to radians
xi=zeta;
beta=(strike+declination);
xd1=X-x0;
xd2=X-(x0 - L/tan(xi)); %%%%% changed - to plus
xd3=X-(x0+w);
xd4=X-(x0 - L/tan(xi)+ w);
%angles and r vectors
phi1=atan2(z,xd1);
r1=sqrt(z^2+xd1.^2);
phi3=atan2(z,xd3); 
r3=sqrt(z^2+xd3.^2);
phi2=atan2(z+L,xd2);
r2=sqrt((z+L)^2+xd2.^2);
phi4=atan2(z+L,xd4);
r4=sqrt((z+L)^2+xd4.^2);

%EQ to model magnetic dike applied to fault location
F1 = 2*k*Fe*sin(xi)*((sin(2*I)*sin(xi)*sin(beta)) - (cos(xi)*((cos(I)^2)*I*(sin(beta))^2 - (sin(I)^2))));
F2 = log((r2.*r3)./(r4.*r1));
F_f = F1.*F2;
F3 = sin(2*I)*cos(zeta)*sin(beta) + sin(zeta)*((cos(I)^2)*(sin(beta)^2) - (sin(I)^2));
F4 = (phi1 - phi2 - phi3 + phi4);
F_l = F3.*F4;
F = F_f+F_l+off;
end

