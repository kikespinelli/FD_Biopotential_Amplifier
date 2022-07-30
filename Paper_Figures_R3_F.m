% An Instrumentation Amplifier INA111 with a gain GINA=51;
% was included at the amplifier differential output
% for measuring with single-ended instruments

close all;
clear all;
GINA=51;



%% GDD
%
%

f1=0.01; f2=10e3;
w1=2*pi*f1; w2=2*pi*f2;

% GDD measurement
% Vin=0.5 mVpp with Vdc=1.5 V
vin=0.5e-3;
f=[0.1 0.2 0.4 0.8 1 2 4 8 10 20 40 80 100 200 400 800 1000];
vo=[10.6 11.1 11.2 11.2 11.2 11.2 11.2 11.2 11.2 11.2 11.1 10.7 10.4 8.72 5.92 3.32 2.72];
Gdd_e=vo/GINA/vin; Gdd_m=mean(Gdd_e(3:10));

% Calculated GDD
% Gn=440
Gn=(1+2*33/2.2)*(1+2*22/3.3); tauL=4.7; tauH=682e-6;

% Frequency response from 0.1 a 1 kHz
% Aproximated expression
Gdd_aprox=zpk([0], [-1/tauH -1/tauL], Gn/tauH);
[mag_a,phase_a,w_a] = bode(Gdd_aprox,{w1,w2});
mag_a=20*log10(squeeze(mag_a));
phase_a=squeeze(phase_a);

% Ganancia de -3dB 
GndB=20*log10(Gdd_m)*ones(length(w_a),1);

% Exact expression
Gdd_x=tf([1 Gn/tauH 0], [1 1/tauH 1/(tauH*tauL)]);
[mag_x,phase_x,w_x] = bode(Gdd_x,{w1,w2});
mag_x=20*log10(squeeze(mag_x));
phase_X=squeeze(phase_x);

% FIGURE
figure(1);

hold off;
% -3db line
p1=semilogx(w_a/(2*pi), GndB-3,'k-.');
hold on;
% Exact curve
p2=semilogx(w_x./(2*pi), mag_x,'Color',[.7 .7 .7],'linewidth',2);
% Measurement
p3=semilogx(f,20*log10(Gdd_e),'ok');
% Approximated curve
p4=semilogx(w_a/(2*pi), mag_a,'k--');

legend([p2,p4,p3,p1], ...
    {'Calculated transfer, exact expression', ...
    'Calculated transfer, approximated expression', ...
    'Measured data points', ...
    '-3 dB line'}, 'Location','SouthWest');
axis([f1 f2 25 55]); 
grid on;
xlabel('Frequency [Hz]');
ylabel('Differential Gain GDD [dB]');

publicationPrint6(gcf,12,10,'Figura1','png',10)


%% Noise spectral density
%
%

% Aproximación de la densidad espectral de ruido del TLC2274
% Está en nV/sqrt(Hz) y tomada de su datasheet
fs=logspace(1,4,1000);

% Experimental Noise data
% El ruido medido está en dos archivos tomados con una ganancia total
% de 22660=51*Gn: 
% Este archivo está armado con dos medidas del Stanford Research SR760: una con
% el rango de 0-50 Hz y otra con 0-1000 Hz. Usé la primera hasta 35 HZ y
% empalmé con la segunda.
load Measured_Noise_TLC2274.txt;
x0=Measured_Noise_TLC2274;

% Calculo Gdd en las frecuencias del espectro
w=x0(:,1)*2*pi;
Gdd_aprox=zpk([0], [-1/tauH -1/tauL], Gn/tauH);
[mag,phase] = bode(Gdd_aprox,w);
mag=squeeze(mag);

% Calculo ruido medido referido a la entrada dividiendo por Gdd
Ei2=(x0(:,2)/GINA).^(2).*mag.^(-2);
ei=Ei2.^(0.5);

% Calculo del ruido total entre 0.5 Hz - 1kHz.
fa=x0(2:533,1);
deltaF=diff(fa); N=531;
Ei2a=Ei2(2:533);
a=0;
for j=1:1:N
    a=a+Ei2a(j)*deltaF(j);
end
RuitoTotal=sqrt(a)
    

% Aproximación de la densidad espectral de ruido del TLC2274
% Está en nV/sqrt(Hz)
fs=w/2/pi;
s=(8+190*fs.^(-0.650))*1e-9; 

% Ruido del resistor R1=3.3k
k=1.38e-23; T=300; Er2=ones(length(s),1)*4*k*T*3.3e3;

% Calculo del ruido total en nV/sqrt(Hz)
st=(2*(s.^2)+Er2).^(0.5);

% Simulación de ruido con TINA
load Noise_TLC2274_TINA.txt; 
load Noise_OPA4344_TINA.txt;

% FIGURE
% Dibujo la Densidad espectral de ruido medida y calculada
% La densidad espectral calculada es válida para f>10Hz que está
% documentada en la datasheet del TLC2274. f=10 Hz corresponde a x0(80)

figure(2);
hold off;
semilogx(w/2/pi, ei*1e9,'k')
hold on;
%semilogx(fs,s*1e9); 
semilogx(w(80:length(w))/2/pi, st(80:length(st))*1e9,'linewidth',2,'Color',[.7 .7 .7]);

semilogx(Noise_TLC2274_TINA(:,1),Noise_TLC2274_TINA(:,2)*1e9,'--k');
semilogx(Noise_OPA4344_TINA(:,1),Noise_OPA4344_TINA(:,2)*1e9,'-.k');

legend('Measured','Estimated (based on datasheet data valid for f>10 Hz)','TLC2274 (simulation)', 'OPA4344 (simulation)')
axis([0.4 1000 0 500]); grid
xlabel('Frequency [Hz]');
ylabel('Noise Spectral Density [nV/\surd Hz]');

publicationPrint6(gcf,12,10,'Figura8','png',10)

%% Transient Response
% Experimental data Pulse_IEC.mat
% Response to a 3 mV 100 ms width pulse. 
% A 3V pulse with a 1000x attenuator was used and a 1.5V battery 
% was added to produce a dc shift
% A periodic 1 Hz pulse train were applied for averaging.
% Agilent Osciloscope in Average mode (32 frames)
% column#1: time, column#2: vin [mV], column#3: vOUT [V]

%%
figure(3);
load Pulse_IEC
t=Pulse_IEC(:,1); vin=Pulse_IEC(:,2); vout=Pulse_IEC(:,3)*1000/440;
hold off;
plot(t,vin-vin(1),'Linewidth',2,'Color',[.7 .7 .7]);
hold on;
plot(t,vout-vout(1),'k'); 
plot([-0.1 1],[-0.1 -0.1],'--k');
axis([-0.1 0.3 -0.2 3.2]);
xlabel('Time [s]'); ylabel('y(t) [mV]');
axis([-0.05 0.20 -0.3 3.3]);
grid on;
publicationPrint6(gcf,12,10,'Figura7','png',10)

%% CMRR Measurements
% An Instrumentation Amplifier INA111 with a gain GINA=101;
% was included at the amplifier differential output
% for measuring with single-ended instruments

% GDC Measurements
% The amplifier was tested using four TLC2274 Integrated Circuits
% Vin=1 Vpp with Vdc=2.5 V
% GDC Trial #1-4

figure(4);
vin=1; GINA=101;
f=[1 2 4 8 10 20 40 80 100 200 400 800 1000];
vo=[428 412 417 398 405 404 396 380 370 315 220 136 116]*1e-3;
Gdc1=vo/(GINA*vin);
vo=[81 82 77 75 76 73 71 72 71 63 51 37 36]*1e-3;
Gdc2=vo/(GINA*vin);
vo=[360 358 362 360 364 363 355 334 317 274 185 101 80]*1e-3;
Gdc3=vo/(GINA*vin);

% GDD measurement
% Vin=0.5 mVpp with Vdc=1.5 V
vin=0.5e-3; GINA=51;
f=[0.1 0.2 0.4 0.8 1 2 4 8 10 20 40 80 100 200 400 800 1000];
vo=[10.6 11.1 11.2 11.2 11.2 11.2 11.2 11.2 11.2 11.2 11.1 10.7 10.4 8.72 5.92 3.32 2.72];
Gdd_e=vo/GINA/vin; Gdd_m=mean(Gdd_e(3:10));

% Common Mode Rejection Ratio CMRR
cmrr1=Gdd_e(5:17)./Gdc1;
cmrr2=Gdd_e(5:17)./Gdc2;
cmrr3=Gdd_e(5:17)./Gdc3;

hold off;
semilogx(f(5:17),20*log10(cmrr1),'k');
hold on;
semilogx(f(5:17),20*log10(cmrr2),'k');
semilogx(f(5:17),20*log10(cmrr3),'k');
semilogx(f(5:17),20*log10(cmrr1),'ok')
semilogx(f(5:17),20*log10(cmrr2),'ok');
semilogx(f(5:17),20*log10(cmrr3),'ok');
axis([1 1000 0 120]); grid;
xlabel('Frequency [Hz]');
ylabel('CMRR [dB]');
publicationPrint6(gcf,12,10,'Figura9','png',10);


%%
% Señales de ECG


figure(5);
load ecg03.txt; x=ecg03;
Ts=1/512; t=0:Ts:Ts*(length(x)-1);

%subtightplot(4,1,2:4,[.07 .01],[.11 .1],[.07 .01]);

subplot(4,1,2:4)
plot(t,-x*1000/440,'k');
xlim([23.2 24.7]);
xlabel('Time [s]'); 
yl1=ylabel('[mV]');
p1 = get(yl1,'pos');
grid on;

subplot(4,1,1)
%subtightplot(4,1,1,[.07 .01],[],[.07 .01]);
load ecg03.txt; x=ecg03;
Ts=1/512; t=0:Ts:Ts*(length(x)-1);
plot(t,-x*1000/440,'k');
xlim([19 33]);
ylim([-1 2]);
yl2=ylabel('[mV]');
p2 = get(yl2,'pos');
set(yl2,'pos',[p2(1)-0.35 p2(2) p2(3)]);
grid on;

publicationPrint6(gcf,12,10,'Figura10','png',10);



