close all
clear
clc

load Referencia_modelo_RCR.mat
p=REF_modelo_RCR(1).P';
q=REF_modelo_RCR(1).Q';
t=REF_modelo_RCR(1).tiempo';

M=length(t);
p(M)=[];
q(M)=[];
t(M)=[];

M=length(t);   % nº de muestras tomadas (nº natural)
tm=1e-3;       % resolución temporal (paso en el tiempo)
fm=1/tm;       % frecuencia de muestreo
to=M*tm;       % periodo de observación (ciclo natural del latido)
fo=1/to;       % resolución espectral (paso en la frecuencia)
f=0:fo:(M-1)*fo;  % vector de frecuencias de la fft (solo la mitad aprox. son positivas)

Q=fft(q);      % componentes de frecuencia de la corriente
P=fft(p-31.7); % componentes de frecuencia de la tensión

Z=P./Q;
plot(f,real(Z),f,imag(Z)); grid on
R1=real(Z(floor(M/2+1)));
R2=real(Z(1))-R1;
C=(-R2-sqrt(R2^2-4*imag(Z(2))^2))/2/imag(Z(2))/2/pi/R2*to;
disp(['R1= ' num2str(R1)])
disp(['R2= ' num2str(R2)])
disp(['C= ' num2str(C)])