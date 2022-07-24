%Range Migration Algorithm from ch 10 of Spotlight
%Synthetic Aperture Radar
%Signal Processing Algorithms, Carrara, Goodman, and Majewski
%***********************************************************************
%a note on formatting, our convention is sif(Xa,t)
% YOU MUST RUN THIS FIRST TO CAL AND BACKGROUND SUBTRACT DATA:
%RMA_FINAL_opendata
%load data
clear all;
load sif;
figcount = 1;
close_as_you_go = 0;
do_all_plots = 0;
set(0,'defaultaxesfontsize',13); %set font size on plots
%so we can see it in the dissertation
% NOTE: the function 'dbv.m' is just dataout = 20*log10(abs(datain));
%***********************************************************************
if do_all_plots == 1
figure(figcount);
S_image = angle(sif); % 取复数的相位角
imagesc(Kr, Xa, S_image);
colormap(gray);
title('Phase Before Along Track FFT');
xlabel('K_r (rad/m)');
ylabel('Synthetic Aperture Position, Xa (m)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'radians','fontsize',13);
print(gcf, '-djpeg100', 'phase_before_along_track_fft.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
end
%**********************************************************************
%along track FFT (in the slow time domain)
%first, symetrically cross range zero pad so that the radar can squint？？？？？？？？？？？？？？？？？？？？
zpad = 1024; %cross range symetrical zero pad  256
szeros = zeros(zpad, size(sif,2));%256*500的0矩阵
for ii = 1:size(sif,2)%1:500
index = (zpad - size(sif,1))/2;
szeros(index+1:(index + size(sif,1)),ii) = sif(:,ii); %symetrical对称  将sif48*500居中,再拓展至256*500
%zero pad
end
sif = szeros;
clear ii index szeros;
S = fftshift(fft(sif, [], 1), 1);  
%S = fftshift(fft(sif, [], 1));
clear sif;
Kx = linspace((-pi/delta_x), (pi/delta_x), (size(S,1)));%因为k=2*pi/lamda k为波数
if do_all_plots == 1
figure(figcount);
S_image = dbv(S);
imagesc(Kr, Kx, S_image, [max(max(S_image))-40,max(max(S_image))]);
colormap(gray);
title('Magnitude After Along Track FFT');
xlabel('K_r (rad/m)');
ylabel('K_x (rad/m)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'dB','fontsize',13);
print(gcf, '-djpeg100', 'mag_after_along_track_fft.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
end

if do_all_plots == 1
figure(figcount);
S_image = angle(S);
imagesc(Kr, Kx, S_image);
colormap(gray);                                                                                                                                       
title('Phase After Along Track FFT');
xlabel('K_r (rad/m)');
ylabel('K_x (rad/m)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'radians','fontsize',13);
print(gcf, '-djpeg100', 'phase_after_along_track_fft.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
end
if do_all_plots == 1
figure(figcount);
S_image = dbv(fftshift(fft(S, [], 2), 2));
imagesc(linspace(-0.5, 0.5, size(S, 2)), Kx, S_image,[max(max(S_image))-40, max(max(S_image))]);
colormap(gray);
title('Magnitude of 2-D FFT of Input Data');
xlabel('R_{relative} (dimensionless)');%无量纲的
ylabel('K_x (rad/m)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'dB','fontsize',13);
print(gcf, '-djpeg100', 'mag_after_2D_fft.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
end
%**********************************************************************
%matched filter
%create the matched filter eq 10.8
aa_max=10;
aa_min=1e3;
for ii = 1:size(S,2) %step thru each time step row to find phi_if
for jj = 1:size(S,1) %step through each cross range in the
%current time step row
%phi_mf(jj,ii) = -Rs*Kr(ii) + Rs*sqrt((Kr(ii))^2 - (Kx(jj))^2);
 aa=sqrt((Kr(ii))^2 - (Kx(jj))^2);%aa即Ky
        if aa>aa_max
            aa_max=aa;
        end
        if aa<aa_min
            aa_min=aa;
        end
phi_mf(jj,ii) = Rs*sqrt((Kr(ii))^2 - (Kx(jj))^2);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Krr(jj,ii) = Kr(ii); %generate 2d Kr for plotting purposes
Kxx(jj,ii) = Kx(jj); %generate 2d Kx for plotting purposes
end
end
smf = exp(j*phi_mf); %%%%%%%%%%%%
%smf = exp(-j*phi_mf); %%%%%%%%%%% THIS IS THE KEY ISSUE !!!!!
%note, we are in the Kx and Kr domain, thus our convention is S_mf(Kx,Kr)

%appsly matched filter to S
S_mf = S.*smf;
%clear smf phi_mf;
if do_all_plots == 1
figure(figcount);
S_image = angle(S);
imagesc(Kr, Kx, S_image);
colormap(gray);
title('Phase After Matched Filter');
xlabel('K_r (rad/m)');
ylabel('K_x (rad/m)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'radians','fontsize',13);
print(gcf, '-djpeg100', 'phase_after_matched_filter.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
end
clear S;
if do_all_plots == 1
figure(figcount);
S_image = dbv(fftshift(fft(S_mf, [], 2), 2));
imagesc(linspace(-0.5, 0.5, size(S_mf, 2)), Kx, S_image,[max(max(S_image))-40, max(max(S_image))]);
colormap(gray);
title('Magnitude of 2-D FFT of Matched Filtered Data');
xlabel('R_{relative} (dimensionless)');
ylabel('K_x (rad/m)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'dB','fontsize',13);
print(gcf, '-djpeg100','mag_after_downrange_fft_of_matched_filtered_data.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
end
%%
%**********************************************************************
%perform the Stolt interpolation
%NOTICE: Must change these parameters!!!!
%Ky_even = linspace(6, 13, 1028); %create evenly spaced Ky for
%interp for book example
%Ky_even = linspace(334.5, 448, 512); %create evenly spaced Ky for
%interp for real data
%Ky_even = linspace(200, 515, 512); %create evenly spaced Ky for
%interp for real data
%FOR DATA ANALYSIS
%kstart = 42.5; %for 1 to 3 ghz
%kstop = 118.5; %for 1 to 3 ghz
kstart =80.3; %for 2 to 4 ghz  85.4
kstop = 103.9; %for 2 to 4 ghz  153.8
%FOR DISSERTATION TO SHOW STOLT WORKING
%kstart = 50;
%kstop = 200;
Ky_even = linspace(kstart, kstop, 512); %create evenly spaced Ky
%for interp for real data
Ky_eeven = linspace(kstart, kstop, zpad); %make this same size as
%kx so we can find downrange
%%%%%%%%%%%%%%%%%clear Ky S_St;

for ii = 1:size(Kx,2) %1:256
Ky(ii,:) = sqrt(Kr.^2 - Kx(ii)^2);%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%S_st(ii,:) = (interp1(Ky(ii,:), S_mf(ii,:), Ky_even)).*H;
S_st(ii,:) = (interp1(Ky(ii,:), S_mf(ii,:), Ky_even));% interp1 差值  256*512
end
S_st(find(isnan(S_st))) = 1E-30; %set all Nan values to 0
%%%%%%%%%%%%%%%%%%%%%%clear S_mf ii Ky;
if do_all_plots == 1
figure(figcount);
S_image = angle(S_st);
imagesc(Ky_even, Kx, S_image);
colormap(gray);
title('Phase After Stolt Interpolation');
xlabel('K_y (rad/m)');
ylabel('K_x (rad/m)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'radians','fontsize',13);
print(gcf, '-djpeg100', 'phase_after_stolt_interpolation.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
end
%*********************************************************************
%perform the inverse FFT's
%new notation: v(x,y), where x is crossrange
%first in the range dimmension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear v Kr Krr Kxx Ky_even;
%v = fftshift(ifft2(S_st,(size(S_st,1)*1),(size(S_st,2)*1)),1);
v = ifft2(S_st,(size(S_st,1)*4),(size(S_st,2)*4));
xx = sqrt(Kx.^2 + Ky_eeven.^2);
%bw = (3E8/(4*pi))*(max(xx)-min(xx));
bw = 3E8*(kstop-kstart)/(4*pi);
max_range = (3E8*size(S_st,2)/(2*bw) - Rs)/.3048;
figure(figcount);
S_image = dbv(v);
imagesc(linspace(-1*(zpad*delta_x/2)/.3048, 1*(zpad*delta_x/2)/.3048, size(v, 1)), linspace(0,-1*max_range, size(v,2)),flipud(rot90(S_image)), [max(max(S_image))-15,max(max(S_image))-0]);
colormap(gray); %MUST DO THIS FOR DISSERTATION FIGS
%colormap('default');
title('Final Image');
ylabel('Downrange (ft)');
xlabel('Crossrange (ft)');
cbar = colorbar;
set(get(cbar, 'Title'), 'String', 'dB','fontsize',13);
print(gcf, '-djpeg100', 'final_image.jpg');
if close_as_you_go == 1
close(figcount);
end
figcount = figcount + 1;
clear cbar close_as_you_go figcount jj do_all_plots;
v = ifft2(S_st); %creat an un-zero padded version of the image
%clear S_st;
save lastimage v max_range zpad delta_x;
%save the set of un-zero padded image data