close all
clear
clc

load Referencia_modelo_RCR.mat
p=REF_modelo_RCR(1).P';
q=REF_modelo_RCR(1).Q';
t=REF_modelo_RCR(1).tiempo';

M=length(t);   % nº de muestras tomadas (nº natural)
tm=1e-3;       % resolución temporal (paso en el tiempo)
fm=1/tm;       % frecuencia de muestreo
to=M*tm;       % periodo de observación (ciclo natural del latido)
fo=1/to;       % resolución espectral (paso en la frecuencia)
f=0:fo:(M-1)*fo;  % vector de frecuencias de la fft (solo la mitad aprox. son positivas)

Q=fft(q);      % componentes de frecuencia de la corriente

R1=0.0256;     % valor de referencia para el modelo
R2=0.4428;     % valor de referencia para el modelo
C=2.1960;      % valor de referencia para el modelo
pout=31.7232;  % valor de referencia para el modelo

Z=R1+R2./(1+1i*2*pi*f(1:floor(M/2+1))*C*R2); % impedancia equivalente de la red RCR
Psim=Q(1:floor(M/2+1)).*Z;  % Construcción de la fft del voltaje: primer tramo
Psim=[Psim conj(Psim(M-floor(M/2):-1:2))];  % Construcción de la fft del voltaje: segundo tramo
psim=real(ifft(Psim))+pout;  % voltaje de salida total
figure
plot(t,p,t,psim)