clear;
clc;
close all;
app=NaN(1);
folder1='C:\Local Matlab Data\2.7GHz-Brian'  %%%%%%%%The folder where you put this Matlab files
cd(folder1)
addpath(folder1)
pause(0.1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Calcultion between ASR-9 and a Base Station at 2.7GHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GMF Info
%%%FAA 841629	U	2715	0		X	Aeronautical Radionavigation	170711	RICHMOND	VA	34	8	ALS	5M00P0N	M1.32	NaN	N/A	RICHMOND	VA	G,ASR-9	G,ASR-9	373020N	0771927W	37.50555556	-77.32416667	373020N	0771927W					M015,IRAC 23320/1,SPS-6285/1	S373	+RADAR PRR SET 27, PRR 1080 AVG,C/W=FAA 841802 +M2770,P/W=FAA 694598 +ATCBI-5 PRR 395 AVG,FA=,EA,ASR,RIC,RIC,RIC,ASR,FTA=EA00992	REQUIRED FOR AIR TRAFFIC CONTROL.	Federal Aviation Administration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Base Station Info
% % % 37 30 49.37N 077 20 30.98W
% % % Their antenna is set to radial 125° with tilt of 3°
% % % Antenna height is 110’
% % % Specific antenna information
% % % Azimuth BW 65°
% % % Elevation BW 7°
% % % Transmit power, as we understand it on the one 100 MHz channel where the edge of the bandwidth is at 2690 MHz,  is 200watts ERP.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Input Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%ITM Inputs
rx_latlon=horzcat(37.50555556,-77.32416667,8);  %%%%%Lat/Lon/Height (meters) [8 meters heights is from the GMF] [Might need to change]
tx_latlon=horzcat(dms2degrees([37 30 49.37]),dms2degrees([-077 20 30.98]),convlength(110,'ft','m'));  %%%%%%%%%%Convert to Decimal Degrees, Convert to Meters: Lat, Lon, Height (meters)
FreqMHz=2700; %%%%%%%Frequency for ITM (MHz)
itm_reliability=50; %%%%% 50%
bs_azimuth=125; % % % Their antenna is set to radial 125° with tilt of 3°
bs_eirp=76; %%%%dBm/100MHz [????]

%%%%%%%%Other Inputs
rx_nf=4.5;  %%%%%%%ASR-9 NF in dB
rx_ant_gain_mb=32; %%%%%%Main Beam gain of ASR-9 in dBi
rx_ant_gain_sl=-22; %%%%%Side Lobe
array_rx_ant_gain=horzcat(rx_ant_gain_mb,rx_ant_gain_sl)
in_ratio=-6; %%%%%I/N Ratio -6dB
radar_threshold=-119.3;

% % % RCVR Noise Power is -113.5 dBm.
% % % Protection Criteria I/N of -6 dB
% % % Maximum Interference power from 5G base station in radar receiver is then -119.3 dBm.


%%%%%%%%%%%%%%%%FDR Inputs
array_rx_if=fliplr(horzcat(0, 0.732, 2, 8))/2; %%%%Frequency MHz (ASR-9) [Need to check these numbers]
array_rx_loss=fliplr(horzcat(0, 3, 20, 60)); %%%%%%%dB Loss (ASR-9)

array_tx_rf=fliplr(horzcat(0, 50, 51, 60)) %%%%Frequency MHz (Base Station) [Half Bandwidth]
array_tx_mask=fliplr(horzcat(0, 0.1, 60, 80)) %%%%%%%dB Loss [Just a placeholder at this point

rx_extrap_loss=-60; %%%%%%%%%RX Extrapolation Slope dB/Decade 60dB (This is generous)
tx_extrap_loss=-60; %%%%%%%%%TX Extrapolation Slope dB/Decade -60dB (This is generous)
rx_freq_mhz=2715;  %%%%%Center Frequency of the ASR-9 from the GMF
tx_freq_mhz=2640; %%%%Center Frequency of the Base Stations
fdr_freq_separation=abs(tx_freq_mhz-rx_freq_mhz)
fdr_calc_mhz=ceil(fdr_freq_separation*1.25)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate ITM Pathloss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
%%%%Initialize
NET.addAssembly(fullfile('C:\USGS', 'SEADLib.dll'));
itmp=ITMAcs.ITMP2P;
TerHandler=int32(1); % 0 for GLOBE, 1 for USGS
TerDirectory='C:\USGS\';

RxLat=rx_latlon(1);
RxLon=rx_latlon(2);
RxHtm=rx_latlon(3);

TxLat=tx_latlon(1);
TxLon=tx_latlon(2);
TxHtm=tx_latlon(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Climate	Derived by using ITU-R.P.617
load('TropoClim.mat','TropoClim')
TropoClim_data=int32(TropoClim);
%%%%%Calculate Radar Radio Climate
tropo_value_rx=get_txt_value_GUI(app,RxLat,RxLon,TropoClim_data); %Gets the Climate of each point
tropo_value=find_tropo_itu617_parfor_GUI(app,horzcat(TxLat,TxLon),TropoClim_data,tropo_value_rx);
RadClim=int32(tropo_value) %%%%%%% 1 Equatorial, 2 Continental Subtorpical, 3 Maritime Tropical, 4 Desert, 5 Continental Temperate, 6 Maritime Over Land, 7 Maritime Over Sea

%%%%%%%%%%%%%%Surface refractivity: Derived by using ITU-R.P.452
load('data_N050.mat','data_N050')
data_N050_data=data_N050;
Refrac=find_refrac_itu452_par_GUI(app,horzcat(TxLat,TxLon),horzcat(RxLat,RxLon),data_N050_data)

%%%%%%%Other ITM inputs
Tpol=1;
Dielectric=25.0;
Conduct=0.02;
ConfPct=50/100;
RelPct=itm_reliability/100;

%%%%%%%%%%%%%%%Do the calculation
[temp_dBloss,propmodeary]=itmp.ITMp2pAryRels(TxHtm,RxHtm,Refrac,Conduct,Dielectric,FreqMHz,RadClim,Tpol,ConfPct,RelPct,TxLat,TxLon,RxLat,RxLon,TerHandler,TerDirectory);
prop_mode=double(propmodeary) %%%%%%%%%% 0 LOS, 4 Single Horizon, 5 Difraction Double Horizon, 8 Double Horizon, 9 Difraction Single Horizon, 6 Troposcatter Single Horizon, 10 Troposcatter Double Horizon, 333 Error
pathloss_dB=double(temp_dBloss)
toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Calculate FDR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[FDR_dB,ED,VD,OTR,DeltaFreq,~,trans_mask_lte]=FDR_ModelII_app(app,fdr_calc_mhz,array_tx_rf,array_rx_if,array_tx_mask,array_rx_loss,tx_extrap_loss,rx_extrap_loss);
toc;

zero_idx=nearestpoint_app(app,0,DeltaFreq);
array_fdr=horzcat(DeltaFreq(zero_idx:end)',FDR_dB(zero_idx:end));

fdr_idx=nearestpoint_app(app,fdr_freq_separation,array_fdr(:,1));
fdr_dB=array_fdr(fdr_idx,:)  %%%%%%Frequency, FDR Loss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FDR Plot
close all;
figure;
hold on;
plot(array_fdr(:,1),array_fdr(:,2),'-b','LineWidth',2,'DisplayName','FDR Loss')
xline(fdr_dB(1),'-g','LineWidth',2,'DisplayName','Delta F [MHz]')
yline(fdr_dB(2),'-r','LineWidth',2,'DisplayName','FDR Loss [dB]')
legend('Location','northwest')
title({strcat('FDR: ASR-9 and Base Station')})
grid on;
xlabel('Frequency Offset [MHz]')
ylabel('FDR [dB]')
filename1=strcat('FDR1.png');
saveas(gcf,char(filename1))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Antennas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Take into consideration the sector/azimuth off-axis gain

%%%%%%%%%%%%%%%%%%%%%%%%%%%ITU Antenna Pattern
deploy_type='macro rural'  %%%%%%%Has 3 degree downtilt
itu_ant_azi=-180:1:180;
temp_theta_range=0;
norm_array_itu_gain_pattern=itu_ant_gain(app,deploy_type,itu_ant_azi,temp_theta_range)-18;  %%%% to normalize to 0dB

%%%%%%%%%Find the azimuth off-axis antenna loss
bs2fed_azimuth=azimuth(tx_latlon(1),tx_latlon(2),rx_latlon(1),rx_latlon(2));  %%%%%Where 0 is North, clockwise.
sector_azi=bs_azimuth;
azi_diff_bs=bs2fed_azimuth-sector_azi;
mod_azi_diff_bs=mod(azi_diff_bs+180,360)-180  %%%%%%%%%%Keep everything within the range of -180 ~ 180
[nn_azi_idx]=nearestpoint_app(app,mod_azi_diff_bs,itu_ant_azi); %%%%%%%Nearest Azimuth Idx
bs_ant_gain=norm_array_itu_gain_pattern(nn_azi_idx)


% % % % %         %%%%%%%%%%Example azimuth calculation with visual
% % % % %         figure;
% % % % %         hold on;
% % % % %         plot(tx_latlon(2),tx_latlon(1),'or')
% % % % %         plot(rx_latlon(2),rx_latlon(1),'sb')
% % % % %         bs2fed_azimuth(1)
% % % % %         grid on;
% % % % %         plot_google_map('maptype','terrain','APIKey','AIzaSyCgnWnM3NMYbWe7N4svoOXE7B2jwIv28F8') %%%Google's API key made by nick.matlab.error@gmail.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Link Budget
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Pr = PT + GT + GR -FDR- Pl – System Loss  (System loss is typical value of 2 dB)
int_margin_dB=bs_eirp+bs_ant_gain+array_rx_ant_gain-fdr_dB(2)-pathloss_dB-2-radar_threshold


%%%int_margin_dB = 11.7466  -42.2534
%%%%%%%%%%%%%%%%%%%%%%%%ASR-9 Main Beam: 11dB over the -6dB I/N
%%%%%%%%%%%%%%%%%%%%%%%%ASR-9 Side Lobe: -42 dB below the -6dB I/N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%This means we need an additional 11dB of loss
   



% % % %%%%%%%%%%%%%%%%%% % % Bob's Notes
% % % % % % ASR-9 3 dB BW 0.564 MHz
% % % % % % ASR-9 NF 3 dB
% % % % % % ASR-9 Antenna Gain 32 dBi.
% % % % % % RCVR Noise Power is -113.5 dBm.
% % % % % % Protection Criteria I/N of -6 dB
% % % % % % Maximum Interference power from 5G base station in radar receiver is then -119.3 dBm.
% % % % % % 
% % % % % % 5G base station Power is a PT of 53.1 dBm.
% % % % % % 5G Bsse station Antenna Gain GT is 23 dBi.
% % % % % % 
% % % % % % Path distance is 1.14 Miles.  That’s close…….
% % % % % % Path loss free space PL is 106 dB.
% % % % % % 
% % % % % % So, we know that:
% % % % % % 
% % % % % % Pr = PT + GT + GR -FDR- Pl – System Loss  (System loss is typical value of 2 dB)
% % % % % % 
% % % % % % Setting Pr = P Imax
% % % % % % 
% % % % % % -119 = 53.10 + 23 + 32 -FDR -PL – Ls
% % % % % % 
% % % % % % Then FDR required to meet the max interference level to meet the protection criteria is:
% % % % % % 
% % % % % % FDR = 119 + 53.1 + 23 + 32 -106 -2
% % % % % % 
% % % % % % FDR = 119 db
% % % % % % 
% % % % % % That’s pretty high. But we can see what it is once we get the emission measurements. Since the 5G base station can’t be filtered, and the radar is not moving. If the base station was moved some it would help the issue or reduce its transmitted power.
% % % % % % Nore the radar is maybe 25 feet or so above ground, and its antenna has a bit up uptilt, the 5G base station antenna is 125 feet high with some down tilt…..so the radar looks right at it when the antenna spins to that azimuth…..Its just all math. Note other places that have this problem have their own geometry of the coupling.

 



















