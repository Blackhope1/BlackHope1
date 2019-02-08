%5. Problema 3 cuerpos: Tierra-Luna-satélite
clc; clear all; close all;
%% DATA
N = 3; %Número de cuerpos  
G = 6.67e-11; %Constante de gravitación universal [N*m^2/kg^2]
m_Earth = 5.972e24; %Masa de la Tierra [kg]
m_Moon = 7.34e22; %Masa de la Luna [kg]
m_sat = 100; %Masa satélite [kg]

mu_Earth = G*m_Earth; %Parámetro gravitacional Tierra [N*m^2/kg]
mu_Moon = G*m_Moon; %Parámetro gravitacional Luna [N*m^2/kg]

%% Órbitas iniciales
%Luna
a_moon = 384.4e6; %Semieje mayor [m]
e_moon = 0.0549; %Excentricidad [-]
i_moon = 23.73*pi/180; %Inclinación respecto al Ecuador [rad]
RAAN_moon = 5*pi/180; %Ascensión recta nodo ascendente [rad]
omega_moon = 180*pi/180; %Argumento del perigeo [rad]
theta_moon = 40*pi/180; %Anomalía verdadera [rad]
%Paso a cartesianos
[r_moon,v_moon] = cartesiano(a_moon,e_moon,i_moon,RAAN_moon,omega_moon,theta_moon,mu_Earth);

%Satélite
%Parámetros orbitales satélite (keplerianos) (iniciales)
rp = 8000e3; ra = 12000e3; 
a_sat = (rp+ra)/2; %Semieje mayor [m]
e_sat = (ra-rp)/(ra+rp); %Excentricidad [-]
i_sat = 50*pi/180; %Inclinación [rad]
RAAN_sat = 20*pi/180; %Ascensión recta nodo ascendente [rad]
omega_sat = 60*pi/180; %Argumento del perigeo [rad]
theta_sat = 10*pi/180; %Anomalía verdadera [rad]
%Paso a cartesianos
[r_sat,v_sat] = cartesiano(a_sat,e_sat,i_sat,RAAN_sat,omega_sat,theta_sat,mu_Earth);

%% Declaración de variables
r = zeros(3,N); %Matriz de vector posición de los N puntos
v = zeros(3,N); %Matriz de vector velocidad de los N puntos

a = input('Elija el tipo de integrador [Euler(1)/RK(2)]: '); %Variable bandera para elegir el tipo de integrador
if a == 1
    m = zeros(N,1); %Vector de masas de los N puntos
    drdt = zeros(3,N); %Matriz derivada de la posición
    dvdt = zeros(3,N); %Matriz derivada de la velocidad

    dij = zeros(3,1); %Distancia entre el punto i y el j
    ai = zeros(3,1); %Aceleración del punto i
end

ntime = 3.8e5; %Instante final de tiempo [s]

if a == 1
    time_step = 100; %Salto temporal [s]
    U = zeros(6,N,ntime/time_step+1); %Matriz de vector de estado de 3 dimensiones
    F = zeros(6,N,ntime/time_step+1); %Matriz derivada del vector de estado de 3 dimensiones
end

%% Inicialización de variables: Caso Tierra-Luna-satélite. Ejes Tierra
r(:,1) = 0.; %Posición Tierra
r(:,2) = r_moon; %Posición inicial Luna en ejes Tierra [m]
r(:,3) = r_sat; %Posición inicial satélite en ejes Tierra [m]
v(:,1) = 0.; %Velocidad Tierra
v(:,2) = v_moon; %Velocidad de rotación Luna alrededor de la Tierra en ejes Tierra [m/s]
v(:,3) = v_sat; %Velocidad de rotación satélite alrededor de la Tierra en ejes Tierra [m/s]

if a == 1
    U(1:3,1,1) = r(:,1); U(1:3,2,1) = r(:,2); U(1:3,3,1) = r(:,3);
    U(4:6,1,1) = v(:,1); U(4:6,2,1) = v(:,2); U(4:6,3,1) = v(:,3);

    m(1) = m_Earth; %Masa Tierra [kg]
    m(2) = m_Moon; %Masa Luna [kg]
    m(3) = m_sat; %Masa satélite [kg]
end

%%  Ley de Fuerzas de Newton
if a == 1
    for k = 1:(ntime/time_step+1)-1
        for i = 1:N
            ai(:) = 0.;
            for j = 1:N
                dij(:) = U(1:3,j,k) - U(1:3,i,k); %Distancia entre los puntos i y j
                if norm(dij)~= 0.
                    ai(:) = ai(:) + G*m(j)*dij(:)/norm(dij)^3; %Aceleración del cuerpo i debido al cuerpo j
                end
            end
            r(:,i) = U(1:3,i,k);
            v(:,i) = U(4:6,i,k);
            drdt(:,i) = v(:,i);
            dvdt(:,i) = ai(:);
            F(1:3,i,k) = drdt(:,i);
            F(4:6,i,k) = dvdt(:,i);

            U(:,i,k+1) = U(:,i,k) + F(:,i,k)*time_step; %Integrador de Euler
        end
    end
else
    [t,U] = ode45(@ecuacion,[0,ntime],[r(:,1);r(:,2);r(:,3);v(:,1);v(:,2);v(:,3)]); %Integrador RK
end

%% Representación gráfica
figure(1) 
if a == 1
    x1 = squeeze(U(1,1,:));
    y1 = squeeze(U(2,1,:));
    z1 = squeeze(U(3,1,:));
    x2 = squeeze(U(1,2,:));
    y2 = squeeze(U(2,2,:));
    z2 = squeeze(U(3,2,:));
    x3 = squeeze(U(1,3,:));
    y3 = squeeze(U(2,3,:));
    z3 = squeeze(U(3,3,:));
    plot3(x1,y1,z1,'b',x2,y2,z2,'r',x3,y3,z3,'g')
else
    plot3(U(:,1),U(:,2),U(:,3),'b',U(:,4),U(:,5),U(:,6),'r',U(:,7),U(:,8),U(:,9),'g')
end
legend('Órbita Tierra','Órbita Luna','Órbita satélite')
xlabel('x[m]')
ylabel('y[m]')
zlabel('z[m]')



function F = ecuacion(t,U)

G = 6.67e-11; %Constante de gravitación universal [N*m^2/kg^2]
m_Earth = 5.972e24; %Masa de la Tierra [kg]
m_Moon = 7.34e22; %Masa de la Luna [kg]
m_sat = 100; %Masa satélite [kg]
m = [m_Earth;m_Moon;m_sat];

F = zeros(18,1);

%Velocidades
F(1:9) = U(10:18);

%Aceleraciones Tierra
F(10) = G*(m(2)/sqrt(((U(4)-U(1))^2+(U(5)-U(2))^2+(U(6)-U(3))^2))^3*(U(4)-U(1)) + m(3)/sqrt(((U(7)-U(1))^2+(U(8)-U(2))^2+(U(9)-U(3))^2))^3*(U(7)-U(1)));
F(11) = G*(m(2)/sqrt(((U(4)-U(1))^2+(U(5)-U(2))^2+(U(6)-U(3))^2))^3*(U(5)-U(2)) + m(3)/sqrt(((U(7)-U(1))^2+(U(8)-U(2))^2+(U(9)-U(3))^2))^3*(U(8)-U(2)));
F(12) = G*(m(2)/sqrt(((U(4)-U(1))^2+(U(5)-U(2))^2+(U(6)-U(3))^2))^3*(U(6)-U(3)) + m(3)/sqrt(((U(7)-U(1))^2+(U(8)-U(2))^2+(U(9)-U(3))^2))^3*(U(9)-U(3)));

%Aceleraciones Luna
F(13) = G*(m(1)/sqrt(((U(1)-U(4))^2+(U(2)-U(5))^2+(U(3)-U(6))^2))^3*(U(1)-U(4)) + m(3)/sqrt(((U(7)-U(4))^2+(U(8)-U(5))^2+(U(9)-U(6))^2))^3*(U(7)-U(4)));
F(14) = G*(m(1)/sqrt(((U(1)-U(4))^2+(U(2)-U(5))^2+(U(3)-U(6))^2))^3*(U(2)-U(5)) + m(3)/sqrt(((U(7)-U(4))^2+(U(8)-U(5))^2+(U(9)-U(6))^2))^3*(U(8)-U(5)));
F(15) = G*(m(1)/sqrt(((U(1)-U(4))^2+(U(2)-U(5))^2+(U(3)-U(6))^2))^3*(U(3)-U(6)) + m(3)/sqrt(((U(7)-U(4))^2+(U(8)-U(5))^2+(U(9)-U(6))^2))^3*(U(9)-U(6)));

%Aceleraciones Sat
F(16) = G*(m(1)/sqrt(((U(1)-U(7))^2+(U(2)-U(8))^2+(U(3)-U(9))^2))^3*(U(1)-U(7)) + m(2)/sqrt(((U(4)-U(7))^2+(U(5)-U(8))^2+(U(6)-U(9))^2))^3*(U(4)-U(7)));
F(17) = G*(m(1)/sqrt(((U(1)-U(7))^2+(U(2)-U(8))^2+(U(3)-U(9))^2))^3*(U(2)-U(8)) + m(2)/sqrt(((U(4)-U(7))^2+(U(5)-U(8))^2+(U(6)-U(9))^2))^3*(U(5)-U(8)));
F(18) = G*(m(1)/sqrt(((U(1)-U(7))^2+(U(2)-U(8))^2+(U(3)-U(9))^2))^3*(U(3)-U(9)) + m(2)/sqrt(((U(4)-U(7))^2+(U(5)-U(8))^2+(U(6)-U(9))^2))^3*(U(6)-U(9)));

end