function [L1,L2,L3,L4,L5,cdm] = PtosLagrange (m1,m2,D)

% Función que calcula el centro de masas del sistema de masas m1 y m2 y...
% obtiene las coordenadas de los puntos de Lagrange referidos a dicho punto.

%Centro de masas respecto del cuerpo más masivo
cdm = (m1*0+m2*D)/(m1+m2);

%Puntos L1, L2 y L3
D1 = cdm; D2 = D-cdm;
%Trabajando en variables adimensionales...
rho = D1/D;
mu1 = 1-rho;
mu2 = rho;

syms x1 x2 x3
ec1 = mu1/((x1+rho)*(x1+rho)) + mu2/(-(x1-1+rho)*(x1-1+rho))-x1 == 0;
ec2 = mu1/((x2+rho)*(x2+rho)) + mu2/((x2-1+rho)*(x2-1+rho))-x2 == 0;
ec3 = mu1/(-(x3+rho)*(x3+rho)) + mu2/(-(x3-1+rho)*(x3-1+rho))-x3 == 0;

S = solve([ec1,ec2,ec3],[x1,x2,x3]);
L1 = [double(S.x1)*D;0;0]; %Primer punto de Lagrange [m]
L2 = [double(S.x2)*D;0;0]; %Segundo punto de Lagrange [m]
L3 = [double(S.x3)*D;0;0]; %Tercer punto de Lagrange [m]

%Puntos L4 y L5
L4 = [(1-2*rho)/2;sqrt(3)/2;0]*D; %Cuarto punto de Lagrange [m]
L5 = [(1-2*rho)/2;-sqrt(3)/2;0]*D; %Quinto punto de Lagrange [m]
end