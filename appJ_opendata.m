%The following MATLAB program was written to open measured SAR data and 
%calibration data. This program calibrates the SAR data and conditions the data to be
%fed into the MATLAB RMA program in appendix C.
%Range Migration Algorithm from ch 10 of Spotlight Synthetic Aperture Radar
%Signal Processing Algorithms, Carrara, Goodman, and Majewski
clear all;
c = 3E8; %(m/s) speed of light
%*********************************************************************
%load IQ converted data here
load rback2 s; %load variable sif %for background subtraction cal data
%*********************************************************************
%perform background subtraction
sif_sub = s;
load rsphere s; %load variable sif %for image data
sif = s-sif_sub; %perform coherent background subtraction
%sif = sif_sub; %image just the background
%sif = s; %image without background subtraction
clear s;
clear sif_sub;
%***********************************************************************
%radar parameters
fc = 2.45E9; %(Hz) center radar frequency (4.069E9 - 1.926E9)/2 + 1.926E9
B = 60E6; %(hz) bandwidth(4.069E9 - 1.926E9)
Tp = 10E-3; %(sec) pulse width
cr = B/Tp; %(Hz/sec) chirp rate
%VERY IMPORTANT, change Rs to distance to cal target
Rs = (37.1)*.3048; %(m) y coordinate to scene center (down range),
%make this value equal to distance to cal target
Xa = 0; %(m) beginning of new aperture length
delta_x = 2*(1/12)*0.3048; %(m) 2 inch antenna spacing
L = delta_x*(size(sif,1)); %(m) aperture length
Xa = linspace(-L/2, L/2, (L/delta_x)); %(m) cross range
%position of radar on aperture L
Za = 0;
Ya = Rs; %THIS IS VERY IMPORTANT, SEE GEOMETRY FIGURE 10.6
t = linspace(0, Tp, size(sif,2)); %(s) fast time, CHECK SAMPLE RATE
Kr = linspace(((4*pi/c)*(fc - B/2)), ((4*pi/c)*(fc + B/2)), (size(t,2)));
%%
%************************************************************************
%callibration
load rcal37pt1_2 s; %load callibration file to standard target
s_cal = s;
load rcalback_2 s; %load background data for cal to standard target
s_cal = s_cal - s; %perform background subtraction
cal = s_cal;
%calculate ideal cal target parameters
%target parameters, 3 targets
at1 = 1; %amplitude of cal target
xt1 = 0;
yt1 = (37.1)*.3048; %(m) distance to cal target
zt1 = 0;
%Rt and Rb for 1 cal target according to equation 10.26
Rb1 = sqrt((Ya - yt1)^2 + (Za - zt1)^2);
xa = 0;
Rt1 = sqrt((xa - xt1).^2 + Rb1^2);
Kr = linspace(((4*pi/c)*(fc - B/2)), ((4*pi/c)*(fc + B/2)), (size(t,2)));
%according to range defined on bottom of page 410
for ii = 1:size(t,2) %step thru each time step to find phi_if
phi_if1(ii) = Kr(ii)*(Rt1 - Rs);
end
cal_theory = at1*exp(-j*phi_if1);
clear phi_if1;
%calculate the calibration factor
cf = cal_theory./(cal);
%apply the cal data
for ii = 1:size(sif,1)
sif(ii,:) = sif(ii,:).*cf; %turn off cal
end
%Save background subtracted and callibrated data
%save sif sif delta_x Rs Kr Xa;
%clear all;
%run IFP
run appC_SBAND_RMA_IFP;