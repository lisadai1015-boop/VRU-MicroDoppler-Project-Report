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

i_plot=0;
for i_frame=1:Ntime
    if i_plot==1
        % ------- plot range-velocity map at differet frames ----------------

        figure(2)
        pcolor(velocity,range,mag2db(abs(squeeze(beamformed(:,:,i_frame)))).');
        hold on;
        plot(velocity(ind_vmax(i_frame)),range(ind_rmax(i_frame)),'ro')
        colormap jet
        shading interp
        %cmax=max(max(abs(squeeze(radar_cube_data(:,:,i_frame))).*window));
        xlabel('Velocity (m/s)')
        ylabel('Range (m)')
        %xlim([-1 1])
        ylim([0 max(range)])
        clim([0 65])
        cid=colorbar('vert');
        ylabel(cid,'Magnitude (dB)')
        title(['Frame #' num2str(i_frame) 'at ' num2str(time(i_frame)) 'seconds'])

        % Save the frame image in the same directory as the pcap file

        if i_movie==1
            if i_frame==1
                exportgraphics(gcf, gifname, Append=false);
            else
                exportgraphics(gcf, gifname, Append=true);
            end
        end
    end
    % ------------ sum up sacttering power near taregt range positions -------
    if i_range_target(i_frame)-nhbin_target>0 && i_range_target(i_frame)+nhbin_target<=Nrange
        V_target(:,i_frame)=sum(abs(beamformed(:,i_range_target(i_frame)-nhbin_target:i_range_target(i_frame)+nhbin_target,i_frame)),2)/(2*nhbin_target+1);
    else
        V_target(:,i_frame)=abs(beamformed(:,ind_rmax(i_frame),i_frame));
    end
    if i_plot==1
        % ----- plot Doppler Velocity at tracked tareget position -----------
        figure(3)
        plot(velocity,V_target(:,i_frame))
        %cmax=max(max(abs(squeeze(radar_cube_data(:,:,i_frame))).*window));
        xlabel('Velocity (m/s)')
        ylabel('Magnitude')
        ylim([0 10])
        title(['Frame #' num2str(i_frame) ' at range of ' num2str(range(i_range_target)) 'm']);
    end
end

% ----------- progagation loss compensation ------------------

% power_factor=(range(ind_rmax)/max(range(ind_rmax(:))));
% power_factor=ones(Nvelocity,1).*power_factor.';
% V_target=V_target.*power_factor; % normalize amplitude to remove range spreading effect

figure(4);subplot(4,1,[1 2])
pcolor(time,velocity,mag2db(V_target));
colormap default
shading flat
cmax3=max(max(mag2db(V_target)));
ylabel('Velocity (m/s)')
xlabel('Time (s)')
ylim([-8 8])
clim([-16.4834 30])
cid=colorbar('northoutside');
ylabel(cid,'Magnitude (dB)')
title('Tracked Target Velocity Responses')

% 

RCS_threshold=std(V_target(:))/10; % nominal 4
buf1=sum(abs(V_target)>RCS_threshold); % sum up target RCS of all velocity

Nsmo=15;
buf_smo=smooth(buf1,Nsmo).';
buf=buf1-buf_smo; %  remove bulk velocity bias 

figure(10);
plot(time,buf1,'g',time,buf_smo,'b',time,buf,'r')
xlabel('Time (s)')
legend('Original','Smoothed','Bias Subtracted')

figure(11);
plot(time,buf,'r')
xlabel('Time (s)')
legend('Original','Smoothed','Bias Subtracted')


% ------- plot target velocity profile for differet frames ----------

% 
% figure(8)
% plot(velocity,std(V_target,[],2),'b','linewidth',2) 
% xlabel('Velocity (m/s)')
% ylabel('Magnitude')
% %legend('STD','Mean','Max')
% title(['Target Velocity STD over ' num2str(max(time)-min(time)) ' Seconds'])

% -----------RCS Normaized Velocity Distribution ----------------

figure(9);subplot(4,1,1)
pcolor(time,velocity,mag2db(V_target));
ax1=gca;
colormap(ax1,'default')
shading flat
ylabel('Velocity (m/s)')
xlabel('Time (s)')
ylim([-8 8])
clim([-10.4628 30])
cid=colorbar('northoutside');
ylabel(cid,'Magnitude (dB)')
title('Tracked Target Velocity')

p=1; % 0.5
threshold=std(V_target(:))*p;
ind=find(V_target(:)>threshold);
V_target_norm=zeros(Nvelocity,Ntime);
V_target_norm(ind)=ones(length(ind),1);

% ---------- calculate velocity spectrum ----------------

buf6=mean((velocity*ones(1,Ntime)).*V_target_norm);
NS=9;
buf6_smo=smooth(buf6,NS).';
buf7=buf6-buf6_smo;

fmax=1/median(diff(time));
NFFT=1024;
md_freq=linspace(0,fmax,NFFT);
win1=tukeywin(Ntime,0.3).';
vd_fft=abs(fft(buf7.*win1,NFFT));
[pks,locs,widths,proms]  = findpeaks(vd_fft(1:NFFT/2),md_freq(1:NFFT/2),'MinPeakProminence',0.2,'SortStr','descend');

figure(9);subplot(4,1,2)
pcolor(time,velocity,V_target_norm);
ax2=gca;
colormap(ax2,'gray')
shading flat
ylabel('Velocity (m/s)')
xlabel('Time (s)')
ylim([-8 8])
title(['Theshold=STD*' num2str(p)])

figure(9);subplot(4,1,3)
buf3=std(V_target_norm,[],2);
plot(velocity,buf3,'b','linewidth',2) 
xlabel('Velocity (m/s)')
ylabel('Magnitude')
ylim([0 0.6])
title('STD')

figure(9);subplot(4,1,4)
if ~isempty(locs)
    xoff=0.04;
    yoff=-0;
    plot(md_freq(1:NFFT/2),vd_fft(1:NFFT/2),'r', locs(1:3),pks(1:3)+0.2,'kv','linewidth',2);
    text(locs(1)+xoff,pks(1)+yoff,[num2str(locs(1)) ' Hz'],'fontsize',12,'color','black');
    text(locs(2)+xoff,pks(2)+yoff,[num2str(locs(2)) ' Hz'],'fontsize',12,'color','black');
    text(locs(3)+xoff,pks(3)+yoff,[num2str(locs(3)) ' Hz'],'fontsize',12,'color','black');
else
    plot(md_freq,vd_fft,'b');
end
ymax=max(vd_fft);
xlabel('Frequency (Hz)')
xlim([0 fmax/2])
ylim([0 ymax*1.5])

figure(10);
plot(time,buf6,'g--',time,buf6_smo,'b-.',time,buf7,'r')
cmax=max(abs(buf6));
xlabel('Velocity (m/s)')
ylabel('Magnitude')
ylim([-cmax*1.5 cmax*1.5])
%title(['Frame #' num2str(i_frame) ' at range of ' num2str(range(i_range_target)) 'm'])
