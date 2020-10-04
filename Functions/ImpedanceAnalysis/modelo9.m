close all
clear
clc

load Sujeto_1D.mat
p=REF(1).P_central';
q=REF(1).Q_central';
t=REF(1).t_central';
pout=REF(1).P_out';

p(length(t))=[]; % elimino la última muestra para que el inicio y el fin no coincidan
q(length(t))=[]; % elimino la última muestra para que el inicio y el fin no coincidan
t(length(t))=[]; % elimino la última muestra para que el inicio y el fin no coincidan

M=length(t);   % nº de muestras tomadas (nº natural)
fm=1000;       % frecuencia de muestreo utilizada
tm=1/fm;       % periodo de muestreo utilizado
to=M/1000;     % periodo de observación de la señal (ciclo de la señal)
fo=1/to;       % resolución espectral (salto de frecuencia tomado por la fft)
f=0:fo:(M-1)*fo;  % vector de frecuencias de la fft (solo la mitad aprox. son positivas)
t=0:tm:(M-1)*tm;  % vector de tiempos

Q=fft(q);      % componentes de frecuencia de la corriente
P=fft(p-pout); % componentes de frecuencia de la tensión en la red

k=0.99;        % coeficiente de energía de la señal
Eq=cumsum(2*((abs(Q(2:floor(M/2+1)))/M).^2)); % vector de energía acumulada con la frecuencia para Q excluyendo la componente continua
Ep=cumsum(2*((abs(P(2:floor(M/2+1)))/M).^2)); % vector de energía acumulada con la frecuencia para P excluyendo la componente continua
Eq=Eq-Eq(floor(M/2+1)-1)*k; % resta del porcentaje de energía total para Q
Ep=Ep-Ep(floor(M/2+1)-1)*k; % resta del porcentaje de energía total para P
Eq=1+sign(Eq); % celdas negativas pasan a 0, positivas a 2
Ep=1+sign(Ep); % celdas negativas pasan a 0, positivas a 2
Eq=find(Eq,1); % encuentro la primera celda no negativa (corresponde a una energía superior al porcentaje)
Ep=find(Ep,1); % encuentro la primera celda no negativa (corresponde a una energía superior al porcentaje)
m=max(Eq,Ep)+1; % determino el número de muestras reales, continua incluida
fmax=(m-1)/to;  % frecuencia máxima a considerar en la señal

Q(m+1:M)=[]; % construyo la nueva fft de Q
P(m+1:M)=[]; % construyo la nueva fft de P
M2=2*m-1;  % el número de muestras ha cambiado
Q=[Q conj(Q((M2-floor(M2/2):-1:2)))]*M2/M; % construyo la nueva fft de Q
P=[P conj(P((M2-floor(M2/2):-1:2)))]*M2/M; % construyo la nueva fft de P  

fm2=M2/to;       % la frecuencia de muestreo ha cambiado
tm2=1/fm2;       % resolución temporal (salto en el tiempo) ha cambiado
f2=0:fo:(M2-1)*fo;  % nuevo vector de frecuencias de la fft (solo la mitad aprox. son positivas)
t2=0:tm2:(M2-1)*tm2;% nuevo vector de tiempos

figure
plot(t2,real(ifft(Q)));hold on; plot(t,q) % aspecto de q(t) con menos muestras
figure
plot(t2,real(ifft(P)));hold on; plot(t,p-pout)% aspecto de p(t) en bornas de la red con menos muestras

figure
subplot(221); plot(f2,abs(Q),'.'); grid on  % Módulo de la fft de la nueva Q
subplot(222); plot(f2,abs(P),'.'); grid on  % Módulo de la fft de la nueva P
subplot(223); plot(f2,angle(Q),'.'); grid on% Fase de la fft de la nueva Q
subplot(224); plot(f2,angle(P),'.'); grid on% Fase de la fft de la nueva P

Z=P./Q;
figure
plot(f2,real(Z),f2,imag(Z)); grid on % real e imaginaria de Z
R1=real(Z(floor(M2/2+1))); % calculo R1
R2=real(Z(1))-R1; % calculo R2
C=(-R2-sqrt(R2^2-4*imag(Z(2))^2))/2/imag(Z(2))/2/pi/R2*to; % calculo C
disp(['R1= ' num2str(R1)])
disp(['R2= ' num2str(R2)])
disp(['C= ' num2str(C)])
disp(['fmax= ' num2str(fmax)])
disp(['energía= ' num2str(100*k)])

Z=R1+R2./(1+1i*2*pi*f2(1:floor(M2/2+1))*C*R2); % impedancia equivalente de la red RCR
Psim=Q(1:floor(M2/2+1)).*Z;  % Construcción de la fft del voltaje: primer tramo
Psim=[Psim conj(Psim(M2-floor(M2/2):-1:2))];  % Construcción de la fft del voltaje: segundo tramo
psim=real(ifft(Psim))+pout;  % voltaje de salida total
figure
plot(t,p,t2,psim);grid on  % comparación entre la presión original y la simulada