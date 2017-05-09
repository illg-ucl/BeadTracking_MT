close all
clear all
clc

snr = [2.1045; 2.011; 1.9582; 1.7928; 1.7466; 1.1575; 0.9883; 0.91739; 0.8247; 0.78969];
mean_gaussian_fitting  = [0.3653;   0.38;   0.41;     0.58;   0.98;   1.3653;  2.37; 2.5443;   2.61;  2.642];
sigma_gaussian_fitting = [0.1894;   0.2041; 0.3685; 0.6651; 0.9969;   1.2453; 1.3276; 1.4707; 1.4622; 1.4968];
mean_mark_fitting      = [0.74;     0.8;    0.82;   0.7842;   0.83;  0.8816;   0.80;   0.87;   0.88; 0.8969;];
sigma_mark_fitting = [0.4297; 0.4480; 0.4580; 0.4302; 0.3887; 0.4303; 0.3649; 0.4025; 0.3872;0.4177];

figure
hold on
plot(snr,mean_gaussian_fitting,'rs-')
plot(snr,mean_mark_fitting,'gs-')
hold off
figure
hold on
plot(snr,sigma_gaussian_fitting,'ro-')
plot(snr,sigma_mark_fitting,'go-')
hold off