clc;
clear all;
close all;
%%MATLAB Codes to Calculate the Optical Power Distribution of LOS Link at Receiving Plane for a Typical Room
theta =70;
% Noise-bandwidth factor %
I2 = 0.562;
% Data rate (Bit per second)
Rb = 115200;
% Ambient light power (Ampere) %
Iamb = 7E-8;
% Photodiode responsivity (A/W )%
R = 0.55;
% Electron charge (C)
q = 1.60E-19;
% Amplifier bandwidth (Hz)%
Ba = 4.5E6;
% Amplifier noise density (Ampere/Hz^0.5)%
Iamf = 5e-12 ;
% semi-angle at half power
ml=-log10(2)/log10(cosd(theta));
%Lambertian order of emission
P_LED=20;
%transmitted optical power by individual LED
nLED=60;
% number of LED array nLED*nLED
P_total=nLED*nLED*P_LED;
%Total transmitted power
Adet=10^-4;
%detector physical area of a PD
Ts=1;
%gain of an optical filter; ignore if no filter is used
index=1.5;
%refractive index of a lens at a PD; ignore if no lens is used
FOV=70;
%FOV of a receiver
G_Con=(index^2)/(sind(FOV).^2);
%gain of an optical concentrator; ignore if no lens is used
%%
lx=10; ly=10; lz=3;
% room dimension in meter
h=2.15;
%led position
x0lc=0;
y0lc=0;
rlc=2;
teta=-pi:0.01:pi;
xlc=rlc*cos(teta)+x0lc;
ylc=rlc*sin(teta)+y0lc;
%hold on
nlc=8;
tet=linspace(-pi,pi,nlc+1);
xilc=rlc*cos(tet)+x0lc;
yilc=rlc*sin(tet)+y0lc;
% %the distance between source and receiver plane
% [XT,YT]=meshgrid([-lx/4 lx/4],[-ly/4 ly/4]);
% % position of LED; it is assumed all LEDs are located at samepoint for
% % faster simulation
% % for one LED simulation located at the central of the room, use XT=0 and YT=0
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=lx*10; Ny=ly*10;
% number of grid in the receiver plane
x=linspace(-lx/2,lx/2,Nx);
y=linspace(-ly/2,ly/2,Ny);
[XR,YR]=meshgrid(x,y);
D1=sqrt((XR-xilc(1,1)).^2+(YR-yilc(1,1)).^2+h^2);
D2=sqrt((XR-xilc(1,2)).^2+(YR-yilc(1,2)).^2+h^2);
D3=sqrt((XR-xilc(1,3)).^2+(YR-yilc(1,3)).^2+h^2);
D4=sqrt((XR-xilc(1,4)).^2+(YR-yilc(1,4)).^2+h^2);
D5=sqrt((XR-xilc(1,5)).^2+(YR-yilc(1,5)).^2+h^2);
D6=sqrt((XR-xilc(1,6)).^2+(YR-yilc(1,6)).^2+h^2);
D7=sqrt((XR-xilc(1,7)).^2+(YR-yilc(1,7)).^2+h^2);
D8=sqrt((XR-xilc(1,8)).^2+(YR-yilc(1,8)).^2+h^2);
% distance vector from source 1
cosphi_A1=h./D1;
cosphi_A2=h./D2;
cosphi_A3=h./D3;
cosphi_A4=h./D4;
cosphi_A5=h./D5;
cosphi_A6=h./D6;
cosphi_A7=h./D7;
cosphi_A8=h./D8;

% angle vector
receiver_angle_1=acosd(cosphi_A1);
receiver_angle_2=acosd(cosphi_A2);
receiver_angle_3=acosd(cosphi_A3);
receiver_angle_4=acosd(cosphi_A4);
receiver_angle_5=acosd(cosphi_A5);
receiver_angle_6=acosd(cosphi_A6);
receiver_angle_7=acosd(cosphi_A7);
receiver_angle_8=acosd(cosphi_A8);


H_A1=(ml+1)*Adet.*cosphi_A1.^(ml+1)./(2*pi.*D1.^2);
H_A2=(ml+1)*Adet.*cosphi_A2.^(ml+1)./(2*pi.*D2.^2);
H_A3=(ml+1)*Adet.*cosphi_A3.^(ml+1)./(2*pi.*D3.^2);
H_A4=(ml+1)*Adet.*cosphi_A4.^(ml+1)./(2*pi.*D4.^2);
H_A5=(ml+1)*Adet.*cosphi_A5.^(ml+1)./(2*pi.*D5.^2);
H_A6=(ml+1)*Adet.*cosphi_A6.^(ml+1)./(2*pi.*D6.^2);
H_A7=(ml+1)*Adet.*cosphi_A7.^(ml+1)./(2*pi.*D7.^2);
H_A8=(ml+1)*Adet.*cosphi_A8.^(ml+1)./(2*pi.*D8.^2);
% channel DC gain for source 1
P_rec_A1=P_total.*H_A1.*Ts.*G_Con;
P_rec_A2=P_total.*H_A2.*Ts.*G_Con;
P_rec_A3=P_total.*H_A3.*Ts.*G_Con;
P_rec_A4=P_total.*H_A4.*Ts.*G_Con;
P_rec_A5=P_total.*H_A5.*Ts.*G_Con;
P_rec_A6=P_total.*H_A6.*Ts.*G_Con;
P_rec_A7=P_total.*H_A7.*Ts.*G_Con;
P_rec_A8=P_total.*H_A8.*Ts.*G_Con;
% received power from source 1;
P_rec_A1(find(abs(receiver_angle_1)>FOV))=0;
P_rec_A2(find(abs(receiver_angle_2)>FOV))=0;
P_rec_A3(find(abs(receiver_angle_3)>FOV))=0;
P_rec_A4(find(abs(receiver_angle_4)>FOV))=0;
P_rec_A5(find(abs(receiver_angle_5)>FOV))=0;
P_rec_A6(find(abs(receiver_angle_6)>FOV))=0;
P_rec_A7(find(abs(receiver_angle_7)>FOV))=0;
P_rec_A8(find(abs(receiver_angle_8)>FOV))=0;
% if the anlge of arrival is greater than FOV, no current is generated at
% the photodiode.
% P_rec_A2=fliplr(P_rec_A1);
% % received power from source 2, due to symmetry no need separate
% % calculations
% P_rec_A3=flipud(P_rec_A1);
% P_rec_A4=fliplr(P_rec_A3);
P_rec_total=P_rec_A1+P_rec_A2+P_rec_A3+P_rec_A4+P_rec_A5+P_rec_A6+P_rec_A7+P_rec_A8;
P_rec_dBm=10*log10(P_rec_total);
Var_power = var(P_rec_dBm);
%%
% 5. For SNR %

Bn = I2 * Rb; % Noise-bandwidth (Sec^-1)%
Pamb = Iamb / R; % Ambient light power (W) %
% Shot-noise variance ( Ampere^2 )%
omega_shot = 2 * q * R * (P_rec_total + Pamb) * Bn; 
% Amplifier noise variance ( Ampere^2 )%
omega_amplifier = Iamf^2 * Ba; 
%Thermal noise variance
omega_thermal = (8*pi*295*112E-8*1E-3*.562*1E6*1.38E-23)+((16*pi^2*1.38E-23*295*1.5*(112E-8)^2*1E-8*.56281E12)/.03);
% Total noise variance ( Ampere^2 )%
omega_total = omega_amplifier + omega_shot+omega_thermal; 
% SNR %
SNR = (( R * P_rec_total )^2)./ omega_total;
SNRdB = 10*log(SNR);
surfc(x,y,P_rec_dBm);
colormap('jet')
% surfc(x,y,capacity);
% contour(x,y,P_rec_dBm);hold on
% mesh(x,y,P_rec_dBm);