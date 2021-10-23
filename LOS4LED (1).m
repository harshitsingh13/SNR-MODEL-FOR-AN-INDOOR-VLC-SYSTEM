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
lx=5; ly=5; lz=3;
% room dimension in meter
h=2.15;
%the distance between source and receiver plane
[XT,YT]=meshgrid([-lx/4 lx/4],[-ly/4 ly/4]);
% position of LED; it is assumed all LEDs are located at samepoint for
% faster simulation
% for one LED simulation located at the central of the room, use XT=0 and YT=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=lx*5; Ny=ly*5;
% number of grid in the receiver plane
x=linspace(-lx/2,lx/2,Nx);
y=linspace(-ly/2,ly/2,Ny);
[XR,YR]=meshgrid(x,y);
D1=sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2);
% distance vector from source 1
cosphi_A1=h./D1;
% angle vector
receiver_angle=acosd(cosphi_A1);
% alternative methods to calculate angle, more accurate if theangle are
% negatives
% nr=[0 0 1];
% RT=[1.25 1.25]; % transmitter location
% for r=1:length(x)
% for c=1:length(y)
%
% angleA12=atan(sqrt((x(r)?1.25).^2+(y(c)?1.25).^2)./h);
% costheta(r,c)=cos(angleA12);
% end
% end
%
%%
% D2=fliplr(D1);
% % due to symmetry
% D3=flipud(D1);
% D4=fliplr(D3);
H_A1=(ml+1)*Adet.*cosphi_A1.^(ml+1)./(2*pi.*D1.^2);
% channel DC gain for source 1
P_rec_A1=P_total.*H_A1.*Ts.*G_Con;
% received power from source 1;
P_rec_A1(find(abs(receiver_angle)>FOV))=0;
% if the anlge of arrival is greater than FOV, no current is generated at
% the photodiode.
P_rec_A2=fliplr(P_rec_A1);
% received power from source 2, due to symmetry no need separate
% calculations
P_rec_A3=flipud(P_rec_A1);
P_rec_A4=fliplr(P_rec_A3);
P_rec_total=P_rec_A1+P_rec_A2+P_rec_A3+P_rec_A4;
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