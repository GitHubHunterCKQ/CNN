%The following MATLAB program was used throughout the dissertation to open 
%simulated SAR data and condition that data to be processed by the RMA MATLAB
%program shown in appendix C.

%Range Migration Algorithm from ch 10 of Spotlight Synthetic
%Aperture Radar
%Signal Processing Algorithms, Carrara, Goodman, and Majewski
%appendixB的代码相当于是  根据appendixA生成的sif 还原出Xa Kr等数据
clear all;
c = 3E8; %(m/s) speed of light
%*********************************************************************
load thruwall s; %load variable sif %for image data
sif = s; %image without background subtraction    sif means signal of intermediate frequency
clear s;
%clear sif_sub;
%***********************************************************************
%radar parameters
fc = 2.45E9; %(Hz) center radar frequency  3G
B = 60E6; %(hz) bandwidth  2E9
Tp = 20E-3; %(sec) pulse width 10ms
cr = B/Tp; %(Hz/sec) chirp rate
%VERY IMPORTANT, change Rs to distance to cal target  将Rs改为到校准目标的距离
Xa = 0; %(m) beginning of new aperture length
delta_x = (2*1/12)*0.3048; %(m) 2 inch antenna spacing!!!!!!  1ft=12inch=0.3048m  相邻两次测量时雷达间距（每次移动距离）
L = delta_x*(size(sif,1)); %(m) aperture length       sif的行数
Xa = linspace(-L/2, L/2, (L/delta_x)); %(m) cross range  %%%相当于还原出appendixA中的Xa

%position of radar on aperture L
Za = 0;
t = linspace(0, Tp, size(sif,2)); %(s) fast time, CHECK SAMPLE RATE 500等分
Kr = linspace(((4*pi/c)*(fc - B/2)), ((4*pi/c)*(fc + B/2)), (size(t,2)));%%%相当于还原出appendixA中的Kr
Rs = 0*.3048; %(m) y coordinate to scene center (down range),
%make this value equal to distance to cal target
Ya = Rs; %THIS IS VERY IMPORTANT, SEE GEOMETRY FIGURE 10.6
%*************************************************e***********************
%Save background subtracted and callibrated data

for ii = 1:size(sif,1)
    sif(ii,:) = sif(ii,:) - mean(sif,1);
end

save sif sif delta_x Rs Kr Xa;
%clear all;
%run IFP

run appC_SBAND_RMA_IFP;