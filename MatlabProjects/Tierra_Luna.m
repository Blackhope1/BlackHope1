%4. Problema 2 cuerpos: Tierra-Luna
clc; clear all; close all; 
%% DATA
N = 2; %Número de cuerpos  
G = 6.67e-11; %Constante de gravitación universal [N*m^2/kg^2]

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

ntime = 1e7; %Instante final de tiempo [s]

if a == 1
    time_step = 100; %Salto temporal [s]
    U = zeros(6,N,ntime/time_step+1); %Matriz de vector de estado de 3 dimensiones
    F = zeros(6,N,ntime/time_step+1); %Matriz derivada del vector de estado de 3 dimensiones
end

%% Inicialización de variables: Caso Tierra-Luna. Ejes Tierra
r(:,1) = 0.; %Posición Tierra
r(1,2) = 384400e3; %Distancia Luna en ejes Tierra [m] según eje x
v(:,1) = 0.; %Velocidad Tierra
v(2,2) = 1e3; %Velocidad de rotación Luna alrededor de la Tierra [m/s]

if a == 1
    U(1:3,1,1) = r(:,1); U(1:3,2,1) = r(:,2);
    U(4:6,1,1) = v(:,1); U(4:6,2,1) = v(:,2);

    m(1) = 5.97e24; %Masa Tierra [kg]
    m(2) = 7.34e22; %Masa Luna [kg]
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
    [t,U] = ode45(@ecuacion,[0,ntime],[r(:,1);r(:,2);v(:,1);v(:,2)]); %Integrador RK
end

%% Representación gráfica
figure(1)
if a == 1
    x1 = squeeze(U(1,1,:));
    y1 = squeeze(U(2,1,:));
    x2 = squeeze(U(1,2,:));
    y2 = squeeze(U(2,2,:));
    plot(x1,y1,'b',x2,y2,'r')
else
    plot(U(:,1),U(:,2),'b',U(:,4),U(:,5),'r')
end
legend('Órbita Tierra','Órbita Luna')
xlabel('x[m]')
ylabel('y[m]')


function F = ecuacion(t,U)

G = 6.67e-11; %Constante de gravitación universal [N*m^2/kg^2]
m_Earth = 5.972e24; %Masa de la Tierra [kg]
m_Moon = 7.34e22; %Masa de la Luna [kg]

m = [m_Earth;m_Moon];

F = zeros(12,1);

%Velocidades
F(1:6) = U(7:12);

%Aceleraciones Tierra
F(7) = G*(m(2)/sqrt(((U(4)-U(1))^2+(U(5)-U(2))^2+(U(6)-U(3))^2))^3*(U(4)-U(1)));
F(8) = G*(m(2)/sqrt(((U(4)-U(1))^2+(U(5)-U(2))^2+(U(6)-U(3))^2))^3*(U(5)-U(2)));
F(9) = G*(m(2)/sqrt(((U(4)-U(1))^2+(U(5)-U(2))^2+(U(6)-U(3))^2))^3*(U(6)-U(3)));

%Aceleraciones Luna
F(10) = G*(m(1)/sqrt(((U(1)-U(4))^2+(U(2)-U(5))^2+(U(3)-U(6))^2))^3*(U(1)-U(4)));
F(11) = G*(m(1)/sqrt(((U(1)-U(4))^2+(U(2)-U(5))^2+(U(3)-U(6))^2))^3*(U(2)-U(5)));
F(12) = G*(m(1)/sqrt(((U(1)-U(4))^2+(U(2)-U(5))^2+(U(3)-U(6))^2))^3*(U(3)-U(6)));

end
