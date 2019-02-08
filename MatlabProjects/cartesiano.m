% Funci�n que devuelve la posici�n y velocidad de un cuerpo que orbita uno
% masivo de constante mu en ejes cartesianos a partir de sus elementos
% orbitales keplerianos

function [rcart,vcart] = cartesiano(a,e,i,RAAN,omega,theta,mu)

% Argin:
%   a: Semieje mayor de la �rbita [m � km]
%   e: Excentricidad [-]
%   i: Inclinaci�n [rad]
%   RAAN: Ascensi�n Recta del Nodo Ascendente [rad]
%   omega: Argumento del perigeo [rad]
%   theta: Anomal�a verdadera [rad]
%   mu: Par�metro gravitacional asociada a la masa en torno a la cual orbita el sat�lite. mu = G*m_cuerpo

%De elementos orbitales a ejes perifocales
%Posici�n
r = (a*(1-e^2))/(1+e*cos(theta)); %Distancia al centro del cuerpo alrededor del cual orbita [[a]]
x = r*cos(theta); y = r*sin(theta); z = 0;
%Velocidad
rdot = sqrt(mu/(a*(1-e^2)))*e*sin(theta); 
vx = -rdot/e; 
if theta == 0 || mod(theta,(2*pi)) == 0
    vy = sqrt(mu/(a*(1-e^2)))*(e+cos(theta));
else
    vy = rdot/(e*sin(theta))*(e+cos(theta));
end
vz = 0;
%Posici�n y velocidad sat�lite en ejes cartesianos
C_Gp = [cos(RAAN)*cos(omega)-sin(RAAN)*cos(i)*sin(omega), -cos(RAAN)*sin(omega)-sin(RAAN)*cos(i)*cos(omega), sin(RAAN)*sin(i); sin(RAAN)*cos(omega)+cos(RAAN)*cos(i)*sin(omega), -sin(RAAN)*sin(omega)+cos(RAAN)*cos(i)*cos(omega), -cos(RAAN)*sin(i); sin(i)*sin(omega), sin(i)*cos(omega), cos(i)]; %Matriz de transformaci�n de ejes perifocales a cartesianos
rcart(1:3,1) = C_Gp*[x;y;z]; vcart (1:3,1) = C_Gp*[vx;vy;vz];
end
