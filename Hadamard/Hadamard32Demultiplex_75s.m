clear all
%% 
AlazarDefs

%File path for saving the data. Change Savefile to 1 to automatically save
%rawdata and reconstructed images.
% FilePath = 'D:\Data\Photothermal_Shared_lab\7-14-21 Gated integration test\'; %Always end with \
% Savefile = 0; %Change to 1 to automatically save all images and rawdata



Image_H_Pixel = 50;
Image_V_Pixel = 50; 

Image_H_frequency = 2/3; %Frequency of the fast-axis trigger in Hz. 2/3 by default (1.5 seconds per fast line, 75 seconds total), change based on pixel dwell time and total imaging time.  
Image_V_frequency = Image_H_frequency/Image_V_Pixel; %Calculates automatically in Hz. Make sure that it equals to slow-axis trigger frequency.
offset = 65500; %The starting point of the first line.Can be tweaked.

%% Specify QCL modulation settings

MidIR_freq = 1; %QCL modulation frequency in kHz (the slowest one when using 2 function generators)
MidIR_period = 1/(MidIR_freq*1000); %QCL modulation period in seconds
                 
%% AlazarCard Acquision Parameters

Sampling_Rate = 1; %0 for 20KSPS 1 for 50KSPS; 2 for 1MSPS; 3 for 10 MSPS; 4 for 25 MSPS; 5 for 50 MSPS;
%Decrease if "buffer overflow" error occurs

% Channel Input Ranges
% 0 for +- 200mV   3 for +- 2V
% 4 for +- 4V      2 for +- 1V - DOESNT WORK WITH THIS CARD
% 1 for +- 400mV   5 for +- 800mv
ChannelARange = 0; 
ChannelBRange = 4; 
ChannalAImpedances = 0; %0 for 50 ohms. 1 for 1 Mohms. Always 0.
ChannalBImpedances = 0; %0 for 50 ohms. 1 for 1 Mohms. Always 0.
External_Trigger = 1; %Don't change unless you know what you are doing
AcqParameters.Timeline = 100; %Determines when the code automatically fails (in seconds). Should be larger than image acquisition time.
TriggerLevel = 20; % In percent of the input range of channel A.

% Switch for sampling rate. Only change if you know what you are doing.
switch Sampling_Rate
    case 0
        ConfigParameters.SamplingRate = SAMPLE_RATE_20KSPS;
        ConfigParameters.GlobalSamplingRate = 20e3;
    case 1
        ConfigParameters.SamplingRate = SAMPLE_RATE_50KSPS;
        ConfigParameters.GlobalSamplingRate = 50e3;
    case 2
        ConfigParameters.SamplingRate = SAMPLE_RATE_1MSPS;
        ConfigParameters.GlobalSamplingRate = 1e6;
    case 3
        ConfigParameters.SamplingRate = SAMPLE_RATE_10MSPS;
        ConfigParameters.GlobalSamplingRate = 10e6;
    case 4
        ConfigParameters.SamplingRate = SAMPLE_RATE_25MSPS;
        ConfigParameters.GlobalSamplingRate = 25e6;
    case 5
        ConfigParameters.SamplingRate = SAMPLE_RATE_50MSPS;
        ConfigParameters.GlobalSamplingRate = 50e6;
    otherwise
        disp('Error! Unsupported sampling rate! Please try number 1-5');
        return;
end

%call this function if the card supports SetBWLimit (see Alazar SDK).
CardSupportBW = 1;

%The number of pre and post trigger samples. The latter should be
%calculated automatically and determines when the acquisition will stop.
%You may change it to a certain value if External_Trigger=0.
PreTriggerSamples = 0; %Usually 0.
PostTriggerSamples = (1.0/Image_V_frequency+1.0/Image_H_frequency)*ConfigParameters.GlobalSamplingRate;

%Don't change this unless you are familiar with Alazar SDK
RecordsPerBuffer = 1; 
BuffersPerAcquisition = 1;

%Number of trigger events to capture. 
NumberofTriggerEvents = 1;

%End of user editable part.
%%
load('\\SimpsonNAS2\Data\SashaRazumtcev\UVF-PTIR\11-8-21 Indomethacin spectrum 75s\hadamard_matrix.mat')
folderpath1 = '\\SimpsonNAS2\Data\SashaRazumtcev\UVF-PTIR\11-8-21 Indomethacin spectrum 75s';
folderlist1 = dir(folderpath1);
filePattern = fullfile(folderpath1, '*.mat');
theFiles = dir(filePattern);
rawimagestack = zeros(50,50,32);
for numberofrawdata = 1:32
    baseFileName = ['Raw_had' num2str(numberofrawdata) '_35us_1kHz_0.85gain_GI.mat'];
    fullFileName = fullfile(theFiles(1).folder, baseFileName);
    load(fullFileName)
    number_of_pixels = (Image_V_Pixel-1)*Image_H_Pixel;
    rawdata1 = rawdata(1:((length(rawdata)/2)-16)); %frequency filtered fluorescence; last 16 points are removed because there is a voltage jump usually
    rawdata2 = rawdata((length(rawdata)/2+1):length(rawdata)); %stage/qcl reference

    %Next part determines the main inputs for data averaging and image
    %reconstruction

    number_of_pixels = Image_V_Pixel*Image_H_Pixel; %total number of pixels in the image
    between_lines = (1/Image_H_frequency)*ConfigParameters.GlobalSamplingRate; %the exact number of datapoints corresponding to one fast axis period
    line_length = 65000; %(1/Image_H_frequency)*ConfigParameters.GlobalSamplingRate*0.9333; %should be a multiple of 50; 0.9333 for 75s/SR1
    pixel_length = 1299; %round(line_length/Image_H_Pixel); %The number+1 should be a multiple of 50

    %This loop divides the rawdata into pixels; all steps are the multiples of
    %the IR period to match the phase on each pixel
    j = 1;
    k = 0;
    pixels_IR = zeros(number_of_pixels,pixel_length+1);
    for i = 1:number_of_pixels
        pixels_IR(i,:) = rawdata1(offset+(j-1)*between_lines+k*pixel_length:offset+(j-1)*between_lines+(k+1)*pixel_length);
        k = k+1;
        if (i+1) - j*Image_H_Pixel > 0
            j = j+1;
            k = 0;
        end
    end

    %Next loop averages the signal on each pixel over IR periods
    samplesperIRperiod = MidIR_period*ConfigParameters.GlobalSamplingRate;
    pixels_avg = zeros(number_of_pixels,samplesperIRperiod);
    for i=1:number_of_pixels
       for j = 1:samplesperIRperiod
           pixel_temp = 0;
           for k = 1:(pixel_length/samplesperIRperiod)
               pixel_temp2 = pixels_IR(i,j+(k-1)*samplesperIRperiod);
               pixel_temp = pixel_temp + pixel_temp2;
           end
           pixels_avg(i,j) = pixel_temp/(pixel_length/samplesperIRperiod);
       end
    end

%     %Calculate modulation depth simply by taking the difference between the
%     %initial fluorescence intensity and intensity during an IR firing event
%     pixels_mod = max(pixels_avg,[],2) - min(pixels_avg,[],2);
% 
%     %Generate and plot image of modulation depth at each pixel
%     Image_modulation = reshape(pixels_mod,[Image_V_Pixel,Image_H_Pixel]);
%     Image_modulation = Image_modulation';
    generate_reference
    phase_initial = find(diff(pixels_ref_avg(1,:)) == max(diff(pixels_ref_avg(1,:))));
    nonlinear_sine_fit_allpix
    rawimagestack(:,:,numberofrawdata) = Image_DLIA;
end
figure;
imagesc(rawimagestack(:,:,1));
colormap(jet);
position = getrect;
position = round(position);
SS2 = squeeze(sum(sum(rawimagestack(position(2):position(2)+position(4),position(1):position(1)+position(3),1:32))));
figure, plot(SS2)
% SS = squeeze(sum(sum(rawimagestack)));
% load('Z:\SashaRazumtcev\UVF-PTIR\10-19-21 Griseofulvin spectrum attempt 1\Hadamard_matrix.mat')
R = SS2'*inv(H_Mat(1:32,1:32));
figure,plot(R)

spectrum = R';
%%
% Image Acq setup
Image_H_Pixel = 200; %200
Image_V_Pixel = 200; %200
Image_H_frequency = 10; % in Hz, 10
Image_V_frequency = Image_H_frequency/Image_V_Pixel; %Image_H_frequency*1/Image_V_Pixel; % in Hz
offset = 6800; %480 for 10/8 and 600 for 5/4 940 for 2.5/2 2300 for 10/8 SR 1
                 
% AlazarCard Setup
Sampling_Rate = 1; %0 for 20KSPS 1 for 1MSPS; 2 for 10MSPS; 3 for 25 MSPS; 4 for 50 MSPS; 5 for 100 MSPS;

% Channel Input Ranges.  
% 0 for +- 200mV        2 for +- 1V 
% 1 for +- 400mV        3 for +- 2V
% 4 for +- 4V           5 for +- 800mv
% 6 for +- 8V at 1MOhms 7 for +- 16V at 1MOhms
ChannelARange = 0; %0 for HWP rotation, 4 for LIA
ChannelBRange = 4; 
ChannalAImpedances = 0; %0 for 50 ohms. 1 for 1 Mohms.
ChannalBImpedances = 0; %0 for 50 ohms. 1 for 1 Mohms.
External_Trigger = 1; %1
AcqParameters.Timeline = 0; %0
% Set up the trigger level in percent.
TriggerLevel = 20; % In percent of the input range of channel A. (-99 to 99)

%Switch for sampling rate input
switch Sampling_Rate
    case 0
        ConfigParameters.SamplingRate = SAMPLE_RATE_20KSPS;
        ConfigParameters.GlobalSamplingRate = 20e3;%20e3;
    case 1
        ConfigParameters.SamplingRate = SAMPLE_RATE_100KSPS;
        ConfigParameters.GlobalSamplingRate = 1e5;
    case 2
        ConfigParameters.SamplingRate = SAMPLE_RATE_10MSPS;
        ConfigParameters.GlobalSamplingRate = 10e6;
    case 3
        ConfigParameters.SamplingRate = SAMPLE_RATE_25MSPS;
        ConfigParameters.GlobalSamplingRate = 25e6;
    case 4
        ConfigParameters.SamplingRate = SAMPLE_RATE_50MSPS;
        ConfigParameters.GlobalSamplingRate = 50e6;
    case 5
        ConfigParameters.SamplingRate = SAMPLE_RATE_100MSPS;
        ConfigParameters.GlobalSamplingRate = 1e8;
    otherwise
        disp('Error! Unsupported sampling rate! Please try number 1-5');
        return;
end

% If the card support SetBWLimit. You have to call this function if the
% card support it...
CardSupportBW = 1;

%Set up the number of pre and post trigger samples.
PreTriggerSamples = 0;
PostTriggerSamples = (1.0/Image_V_frequency+1.0/Image_H_frequency)*ConfigParameters.GlobalSamplingRate; %10e6;

%Advance number of trigger events settings.
%In order to make alazar card happy, you need to make sure you have enough
%samplings in one buffer. A reference point for the total number of
%samplings in a buffer is 2048*100 = 204800. 
RecordsPerBuffer = 1; %% not sure what this means
BuffersPerAcquisition = 1; 

%Number of trigger events you want to caputure.!!!Auto Setting, not working
%at this stage. Please use advanced settings.
NumberofTriggerEvents = 1;

%End of user editable part.
folderpath1 = '\\SimpsonNAS2\Data\SashaRazumtcev\UVF-PTIR\10-21-21 Griseofulvin spectrum attempt 2\Fluorescence images';
folderlist1 = dir(folderpath1);
filePattern = fullfile(folderpath1, '*.mat');
theFiles = dir(filePattern);
fluimagestack = zeros(200,200,length(theFiles));
for numberofrawdata = 1:length(theFiles)+1
    baseFileName = ['fluor_had1' num2str(numberofrawdata) '.mat'];
%     if numberofrawdata == 22
%         baseFileName = ['fluor__had1' num2str(23) '.mat'];
%     end
    fullFileName = fullfile(theFiles(1).folder, baseFileName);
    load(fullFileName)
    rawdata2 = rawdata((length(rawdata)/2+1):length(rawdata));
    rawdata1 = rawdata(1:length(rawdata)/2);

    % rawdata1 = rawdata1 - 32700;
    % for i = 1:length(rawdata1)
    %         rawdata1(i) = abs(rawdata1(i));
    % end

    samplesperline = (1.0/Image_H_frequency)*ConfigParameters.GlobalSamplingRate;
    data = rawdata1(offset+1:offset+(1.0/Image_V_frequency)*ConfigParameters.GlobalSamplingRate);
    data = reshape(data,[samplesperline,Image_V_Pixel]);

    data2 = rawdata2(offset+1:offset+(1.0/Image_V_frequency)*ConfigParameters.GlobalSamplingRate);
    data2 = reshape(data2,[samplesperline,Image_V_Pixel]);

    index = ones(Image_H_Pixel+1,1);
    for i = 1:Image_H_Pixel
        x = acos(1 - 2*i/Image_H_Pixel);
        index(i+1) = fix(x/pi*(samplesperline/2));
    end 

    RawImage = ones(Image_V_Pixel,Image_H_Pixel);
    RawImage2 = ones(Image_V_Pixel,Image_H_Pixel);
    for i = 1:Image_V_Pixel
        line_data = (data(1:samplesperline/2,i));%+flip(i,data(samplesperline/2+1:end)))/2;
        for j = 1:Image_H_Pixel                      
            RawImage(i,j) = mean(line_data(index(j):index(j+1)));
        end
    end

    fluimagestack(:,:,numberofrawdata) = RawImage;
    
end
RI = squeeze(sum(sum(fluimagestack)));
R2 = (SS2./RI)'*inv(H_Mat(1:31,1:31));
figure,plot(R2)