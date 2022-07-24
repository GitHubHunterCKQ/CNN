%The following MATLAB program was written to simulate SAR data of three 
%different point scatterers at various locations in a target scene.

%Range Migration Algorithm from ch 10 of Spotlight Synthetic
%Aperture Radar
%Signal Processing Algorithms, Carrara, Goodman, and Majewski
clear all;
c = 3E8; %(m/s) speed of light
%**********************************************************************
%radar parameters
fc = 2.45E9; %(Hz) center radar frequency  3
B = 60E6; %(hz) bandwidth 2G
Rs = 0;%10*.3048; %(m) y coordinate to scene center (down range)  Rs为参考距离
Xa = 0; %(m) beginning of new aperture length
L = 8*.3048; %(m) aperture length  8ft  1ft=0.3084m
Xa = linspace(-L/2, L/2, 48); %(m) cross range position of radar on aperture L   %选择48个位置
Za = 0*.3048;
Ya = 0; %THIS IS VERY IMPORTANT, SEE GEOMETRY FIGURE 10.6
fsteps = 5000;%每个位置采样点数
%**********************************************************************
%create SAR if data according to eq 10.4 and 10.5 (mocomp
%to a line) ignoring RVP term
%target parameters, 3 targets
at1 = 1;
xt1 = 5*.3048;
yt1 = -10*.3048;
zt1 = 0;
at2 = 1;
xt2 = 0*.3048;
yt2 = -25*.3048;
zt2 = 0;
at3 = 1;
xt3 = -5*.3048;
yt3 = -10*.3048;
zt3 = 0;
%Rt and Rb for 3 targets according to equation 10.26
Rb1 = sqrt((Ya - yt1)^2 + (Za - zt1)^2);
Rt1 = sqrt((Xa - xt1).^2 + Rb1^2);%  三维坐标系中目标到雷达的直线距离
%Rt1 = sqrt(yt1^2 + (Xac-xt1).^2);
Rb2 = sqrt((Ya - yt2)^2 + (Za - zt2)^2);
Rt2 = sqrt((Xa - xt2).^2 + Rb2^2);
%Rt2 = sqrt(yt2^2 + (Xac-xt2).^2);
Rb3 = sqrt((Ya - yt3)^2 + (Za - zt3)^2);
Rt3 = sqrt((Xa - xt3).^2 + Rb3^2);
%Rt3 = sqrt(yt3^2 + (Xac-xt3).^2);
Kr = linspace(((4*pi/c)*(fc - B/2)), ((4*pi/c)*(fc + B/2)), fsteps);%相当于每个位置5000个数据
%according to range defined on bottom of page 10
%s(xn,w(t))=at*e^(-j*Kr*(Rt-Rs)) Kr=2*w(t)/c     !!!!!!!!!!!!!!!!! delta_phi=4*pi*delta_d/lamda !!!!!!!
%%
for ii = 1:fsteps %step thru each time step to find phi_if  500
for jj = 1:size(Xa,2) %step thru each azimuth step  size(Xa,2)=Xa的列数  48
phi_if1(jj,ii) = Kr(ii)*(Rt1(jj) );
phi_if2(jj,ii) = Kr(ii)*(Rt2(jj) );
phi_if3(jj,ii) = Kr(ii)*(Rt3(jj) );
end
end
sif1 = at1*exp(-j*phi_if1);
sif2 = at2*exp(-j*phi_if2);
sif3 = at3*exp(-j*phi_if3);
sif = sif1+sif2+sif3; %superimpose all three targets
% clear sif1;
% clear sif2;
% clear sif3;
% clear phi_if1;
% clear phi_if2;
% clear phi_if3;
s = sif;
save thruwall s;
%view simulated range history after range compression (figure 10.7) OK
%for ii = 1:size(sif,1);
%  sview(ii,:) = fftshift(fft(sif(ii,:)));
%end
%***********************************************************************
%a note on formatting, our convention is sif(Xa,t)




%**************************************************************
%plot the real value data for the dissertation
set(0,'defaultaxesfontsize',13);
imagesc(Kr*c/(4*pi*1E9), Xa/.3048, real(sif));%Kr*c/(4*pi*1E9)=fc+-B
colormap(gray);
ylabel('x position (ft)');
xlabel('recieved chirp frequency (GHz)');
title('real values of the single point scatter SAR data matrix');
colorbar;
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'V/m','fontsize',13);
