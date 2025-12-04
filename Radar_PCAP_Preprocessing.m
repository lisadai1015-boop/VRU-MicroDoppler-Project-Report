clear all
[fileName, filePath] = uigetfile('*.pcap', 'Select the pcap file');
if fileName == 0
    error('No file selected. Exiting script.');
end
Name = fullfile(filePath, fileName);

pcapReaderObj = pcapReader(Name);
decodedPackets = readAll(pcapReaderObj);
clear pcapReaderObj
NP=size(decodedPackets,2); % total number of packets
total_damaged_frames = 0; %total damaged/incomplete frames
Data = [];  % Array to store the radar data prior reshaping
Buf = []; % Array to store temporary radar data
frameCounter = 0;  % Frame Counter to indicate how many full frames we get
frameDataCounter = 0; % Frame Data Counter / Datapoints per frame
temp = []; %Test array to store the message length after the headers
flagCounter=[]; %Flag Counter Array
time=[];
cube.TimeStamp=[];

%% Check to make sure we start collecting the data when frame flag is 1
flag_1=0;
flag_1st_Frame=0;
for m=1:NP
    bytebuf=dec2hex(decodedPackets(m).Packet.eth.Payload);
    buf2=transpose(bytebuf);
    buf2=buf2(:).';
    ind=strfind(buf2,'0B1103E8C355');
    offset=20;
    if ~isempty(ind)
        buf3=buf2(ind(1)+offset:end);
        buf3=transpose(reshape(buf3,2,[]));
        cube.Bytes=uint8(hex2dec(buf3));
        cube = SMS_Transport_Protocol(cube);    % ************* Read the SMS header
        cube = Debug_Port_Header(cube);         % ************* Read the Debug header
        flagCounter = [flagCounter;cube.FrameFlags]; %///////////////

        if cube.FrameFlags == 1  %FBegining of Framee
            Buf = []; % Empty buf array to store the next frame data
            frameDataCounter = 0;  % Clear the data counter for the next fram
            if flag_1st_Frame==0
                fprintf('First frame detected \r');
                flag_1st_Frame=1; % mark begining of the first frame
            end
            fprintf(['Reading frame #' num2str(frameCounter+1) ' started --- ']);

            cube = Generic_Port_Header(cube); % ************** Read the generic port header
            cube = Static_Port_Header(cube);  % ************** Read the static port header

            if flag_1==0  % get radar parameter once
                sizeOfData = uint32(cube.doppler_bins) * uint32(cube.rx_channels) * uint32(cube.range_gates) * uint32(cube.chirp_types); % Data size of the frame
                doppler_bins = double(cube.doppler_bins);
                rx_channels = double(cube.rx_channels);
                range_gates = double(cube.range_gates);
                sequence = double(cube.chirp_types);
                flag_1=1;
            end
            [Buf, frameDataCounter] = Radar_Data_Decode(cube, Buf, 87, frameDataCounter); % Read radar cube data from the starting index 87

        elseif cube.FrameFlags == 0 && flag_1st_Frame==1 %Frame flag is 0 (Continue of Frame)
            [Buf, frameDataCounter] = Radar_Data_Decode(cube, Buf, 23, frameDataCounter); % Read radar cube data from the starting index 23

        elseif cube.FrameFlags == 2 && flag_1st_Frame==1 %Frame flag is 2 (End of Frame)
            [Buf, frameDataCounter] = Radar_Data_Decode(cube, Buf, 23, frameDataCounter); % Read radar cube data from the starting index 23

            fprintf(['Reading frame #' num2str(frameCounter+1) ' completed \r'])

            % Check to see if the data size matches what we expect
            if ~(frameDataCounter == sizeOfData)
                fprintf("%.2f / %.2f \n",frameDataCounter,sizeOfData);
                fprintf("Incorrect frame data size \r");
                total_damaged_frames = total_damaged_frames + 1;
                Buf = []; % Empty buf array to store the next frame data
                frameDataCounter = 0;  % Clear the data counter for the next frame
                %break
            else
                Data = [Data, Buf];  % Append the frame data into the radar data array
                time=[time; cube.TimeStamp(end)/1e6];
                frameCounter = frameCounter + 1;  % Increment frame counter
                Buf = []; % Empty buf array to store the next frame data
                frameDataCounter = 0;  % Clear the data counter for the next fram
            end % end of data frame size checking
        end% end of FrameFlags checking
    end % End '7E' pattern checking
end % End NP packet loop
fprintf("%f",total_damaged_frames)

%% Data reshaping and process

% 5D array that saves the radar cube data accordingly

% //////// comment out next line if data are loaded from .mat //////////
radar_cube_data = reshape(Data, doppler_bins, rx_channels, range_gates, sequence, frameCounter);
i_SQ=1;     % keep SQ A data only
radar_cube_data = squeeze(radar_cube_data(:,:,:,i_SQ,:));
radar_cube_data=fftshift(radar_cube_data,1);  % move zero velocity bin to middle of array
radar_cube_data=double(radar_cube_data);
% ------------ generate range vector -----------------------

rfftRangePerBinRes=0.332453;    % in meters
vfftSpeedPerBinRes=0.129992;  % in m/s
rfftBinPerSpeedRes= -0.0261659;  % unitless
range=transpose(0:rfftRangePerBinRes:rfftRangePerBinRes*(range_gates-1));

% --------- generate doppler velocity vector -----------------------

velocity=transpose(0:vfftSpeedPerBinRes:vfftSpeedPerBinRes*(doppler_bins-1));
velocity=velocity-(max(velocity)+min(velocity))/2;
time=time-time(1);

% //////// keep only SQ 1 and range&velocity regions of interest ////////

r_min=0;    % specify minimum distance in meteres to keep
r_max=30;   % specify maximum distance in meters to keep
v_max=8;    % specify maximum doppler velopcty in m/s to keep

ind_r=find(range>=r_min & range<=r_max);
ind_v=find(abs(velocity)<=v_max);
range=range(ind_r);
velocity=velocity(ind_v);
radar_cube_data= squeeze(radar_cube_data(ind_v,:,ind_r,:,:)); % apply gating

% /////// Preview RX1 R&V farmes for data quality checking /////////

i_Rx=1;  % keep Rx1 only
radar_cube_data_1CH=double(squeeze(radar_cube_data(:,i_Rx,:,:)));

figure(2)
for i_frame=1:frameCounter
    pcolor(velocity,range,mag2db(abs(squeeze(radar_cube_data_1CH(:,:,i_frame)))).');
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
    % fig=gcf;
    % F(i_frame)=getframe(fig);
    % exportgraphics(gcf, gifFile, Append=true);
    pause(0.1);
end

% ------------- plot time profile -----------

figure(1)
plot(time(2:end),diff(time))
ylabel('Frame Rate (s)')
xlabel('Time (s)')
% ------------- same .mat file ----------------

ind=findstr(Name,'.pcap');
matname=[Name(1:ind-1) '_preprocessed.mat'];
save(matname,'velocity','range','time','frameCounter','radar_cube_data'); 