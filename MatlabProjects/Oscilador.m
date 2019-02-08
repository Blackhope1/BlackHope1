%1. Oscilador arm�nico
clc; clear all; close all;

%% Elegir integrador
b = input('Seleccione c�digo-n�mero del integrador temporal indicado entre par�ntesis [Euler(1)/EInverso(2)/RK(3)]: ');

%% Declaraci�n de variables
x = 1.; %Posici�n inicial
v = 0.; %Velocidad inicial

ntime = 50; %Instante final de tiempo [s]
time_step = 0.01; %Salto temporal [s]
U = zeros(2,ntime/time_step+1); %Matriz de vector de estado
if b == 1 || b == 2
    F = zeros(2,ntime/time_step+1); %Matriz derivada del vector de estado
end

%% Inicializaci�n de variables
U(1,1) = x; U(2,1) = v;

%% Integrador temporal
if b == 1
    for k = 1:1:(ntime/time_step+1)-1
        F(1,k) = U(2,k);
        F(2,k) = -U(1,k);

        U(:,k+1) = U(:,k) + F(:,k)*time_step; %Integrador de Euler
    end
elseif b == 2
    for k = 1:1:(ntime/time_step+1)-1
        F(1,k) = U(2,k);
        F(2,k) = -U(1,k);

        U(:,k+1) = (eye(2)-time_step*[0, 1;-1,0])\U(:,k); %Integrador de Euler inverso
    end
elseif b == 3
    [t,U] = ode45(@ecuacion,[0,ntime],[x;v]); %Integrador RK
end

%% Soluci�n anal�tica
if b == 1 || b == 2
    x_ana = cos(0:time_step:ntime);
    error = abs(x_ana - U(1,:));
elseif b == 3
    x_ana = cos(t);
    error = abs(x_ana - U(:,1));
end


%% Representaci�n gr�fica
figure(1) %Posici�n num�rica
if b == 1 || b == 2
    plot(0:time_step:ntime,U(1,:),'-',0:time_step:ntime,x_ana,'--')
elseif b == 3
    plot(t,U(:,1),'-',t,x_ana,'--')
end
xlabel('Time');
ylabel('X');
title('Solution to Harmonic Equation');
legend('Numerical','Analytical')

figure(2) %Error en diferencias entre la soluci�n num�rica y la anal�tica
if b == 1 || b == 2
    plot(0:time_step:ntime,error)
elseif b == 3
    plot(t,error)
end
xlabel('Time');
ylabel('Error');
title('Solution to Harmonic Equation Error');


function DerivadaEstado = ecuacion(t,U)
%Ecuaci�n del movimiento oscilador arm�nico
DerivadaEstado = [U(2);-U(1)];
end