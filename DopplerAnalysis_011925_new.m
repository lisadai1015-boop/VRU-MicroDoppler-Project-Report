clear all
[fileName, filePath] = uigetfile('*.mat', 'Select the mat file');
if fileName == 0
    error('No file selected. Exiting script.');
end
Name = fullfile(filePath, fileName);
ind=findstr(Name,'.mat');
gifname=[Name(1:ind) 'gif'];

% load data
load(Name); 

[Nvelocity, Nrx, Nrange, Ntime]=size(radar_cube_data);
radar_cube_data=permute(radar_cube_data,[1 3 2 4]);

% -------- calculate max accross different frame as background ----------------

radar_std=squeeze(std(radar_cube_data,[],4));       % max over frames
radar_std=squeeze(std(radar_std,[],3));             % max over Rx channel

Vzero_blanking=ones(Nvelocity,Nrange,Nrx);
v_blanking=0.005;
dr= range(2)-range(1);
dv=abs(velocity(2)-velocity(1));
NH_blanking=ceil(v_blanking/dr);
index_blanking=(round(Nvelocity/2)-NH_blanking+1:round(Nvelocity/2)+NH_blanking+1).';

for m_rx=1:Nrx
    Vzero_blanking(index_blanking,:,m_rx)=zeros(length(index_blanking),1)*ones(1,Nrange);
end

w_range=1;              % max extent in meters of moving parts
Nw=ceil(w_range/dr);

irx_ref=1;

beamformed=zeros(Nvelocity,Nrange,Ntime);
ind_vmax=zeros(Ntime,1);
ind_rmax=zeros(Ntime,1);

% --------------------
%Vzero_blanking=1;
for m_frame=1:Ntime
    % ------------- find target position ----------------
    rd_buf=squeeze(radar_cube_data(:,:,:,m_frame)).*Vzero_blanking;
    rd_buf=reshape(rd_buf,Nvelocity*Nrange,Nrx);
    [buf,indmax]=max(abs(rd_buf(:,irx_ref))); 

    ind_vmax(m_frame)=mod(indmax,Nvelocity);
    ind_rmax(m_frame)=floor(indmax/Nvelocity)+1;
    % ---------- Combine Rx Channels to Beamforming to broadside ------------

    A=squeeze(rd_buf(indmax,:)).';
    Wi=(rand(Nrx,1)*2-1)*pi; % initial phpase weight
    window=ones(Nrx,1);  % initial amplitude weight
    %options = optimset('Display','off','PlotFcns',@optimplotfval,'TolFun',1e-7,'TolX',1e-7);
    options = optimset('Display','off','PlotFcns',[],'TolFun',1e-7,'TolX',1e-7);
    W = fminsearch(@PhaseWeightedSum,Wi,options,[A window]); % find optimal phase weightings for maximizing array gain for each direction
    rd_buf=reshape(rd_buf,Nvelocity,Nrange,Nrx);
    for m_Rx=1:Nrx
        beamformed(:,:,m_frame)=beamformed(:,:,m_frame)+squeeze(rd_buf(:,:,m_Rx))*exp(1i*W(m_Rx))/Nrx;
    end
end


%% ---- power max over time (Untracked) ---------

RadardBMax=max(mag2db(abs(beamformed)),[],3);

figure(1);subplot(131)
pcolor(velocity,range,RadardBMax.');
colormap default
shading flat
cmax2=max(RadardBMax(:));
xlabel('Velocity (m/s)')
ylabel('Range (m)')
clim([cmax2-30 cmax2])
cid=colorbar('vert');
ylabel(cid,'Magnitude (dB)')
title('Max oevr Time (Untracked)')

%% -------------- power sum over time (Untracked) ---------

RadarPowerSum=sum(abs(beamformed).^2,3);

figure(1);subplot(132)
pcolor(velocity,range,pow2db(RadarPowerSum).');
colormap default
shading flat
cmax2=max(pow2db(RadarPowerSum(:)));
xlabel('Velocity (m/s)')
ylabel('Range (m)')
clim([cmax2-30 cmax2])
cid=colorbar('vert');
ylabel(cid,'Magnitude (dB)')
title('Power Sum oevr Time (Untracked)')

%% -------------- STD over time (Untracked) ---------

RadarSTD=std(beamformed,[],3);
figure(1);subplot(133)
pcolor(velocity,range,mag2db(RadarSTD).');
colormap default
shading flat
cmax2=max(mag2db(RadarSTD(:)));
xlabel('Velocity (m/s)')
ylabel('Range (m)')
clim([cmax2-30 cmax2])
cid=colorbar('vert');
ylabel(cid,'Magnitude (dB)')
title('STD oevr Time (Untracked)')


% ---- plot beamformed tracked range-doppler response accross all frames ---------
i_movie=0;
V_target=zeros(Nvelocity,Ntime);
dr_target=2;
nhbin_target=ceil(dr_target/dr/2);

i_move=input('(1) Moving target (0)Stationary Target: ')
if i_move==0
    i_range_target=round(mean(ind_rmax))*ones(length(ind_rmax),1);    % uncomment for fixed position
else
    i_range_target=ind_rmax;
end

