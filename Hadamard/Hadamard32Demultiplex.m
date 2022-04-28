clear all
%% 
AlazarDefs

%File path for saving the data. Change Savefile to 1 to automatically save
%rawdata and reconstructed images.
FilePath = 'D:\Data\Photothermal_Shared_lab\7-14-21 Gated integration test\'; %Always end with \
Savefile = 0; %Change to 1 to automatically save all images and rawdata

IRChannelNumber = 97; %Indicate the used QCL channel for the file name. (Only needed when Savefile = 1)


% Specify image acquisition parameters

%Number of pixels on each axis. 100x100 by default to increase spectral
%accuracy. Use 200x200 for faster high-resolution imaging if SNR permits
Image_H_Pixel = 100;
Image_V_Pixel = 100; 

Image_H_frequency = 1/3; %Frequency of the fast-axis trigger in Hz. 1/3 by default (3 seconds per fast line, 300 seconds total), change based on pixel dwell time and total imaging time.  
Image_V_frequency = Image_H_frequency/Image_V_Pixel; %Calculates automatically in Hz. Make sure that it equals to slow-axis trigger frequency.
offset = 480; %Artifact from the older versions. Not used, but keep this line to avoid errors.

% Specify QCL modulation settings

MidIR_freq = 1; %QCL modulation frequency in kHz (the slowest one when using 2 function generators)
MidIR_period = 1/(MidIR_freq*1000); %In seconds
                 
% AlazarCard Acquision Parameters

Sampling_Rate = 1; %0 for 20KSPS 1 for 50KSPS; 2 for 1MSPS; 3 for 10 MSPS; 4 for 25 MSPS; 5 for 50 MSPS;
%Decrease if "buffer overflow" error occurs

% 0 for +- 200mV   3 fo
% Channel Input Rangesr +- 2V
% 4 for +- 4V           2 for +- 1V - DOESNT WORK WITH THIS CARD
% 1 for +- 400mV             5 for +- 800mv
ChannelARange = 0; %Usually 4 (largest input range)
ChannelBRange = 4; %Almost always 4 for Lock-In detection.
ChannalAImpedances = 0; %0 for 50 ohms. 1 for 1 Mohms. Always 0.
ChannalBImpedances = 0; %0 for 50 ohms. 1 for 1 Mohms. Always 0.
External_Trigger = 1; %
AcqParameters.Timeline = 360; %In seconds. Should be larger than image acquisition time. Use 360 for 300s acquisition.
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

% If the card support SetBWLimit, you have to call this function.
CardSupportBW = 1;

%The number of pre and post trigger samples. The latter should be
%calculated automatically and determines when the acquisition will stop.
%You may change it to a certain value if External_Trigger=0.
PreTriggerSamples = 0; %Usually 0.
PostTriggerSamples = (1.0/Image_V_frequency+1.0/Image_H_frequency)*ConfigParameters.GlobalSamplingRate; %10e6;

%Advanced number of trigger events settings.
%In order to make the alazar card happy, you need to make sure you have enough
%samplings in one buffer. A reference point for the total number of
%samplings in a buffer is 2048*100 = 204800. 
RecordsPerBuffer = 1; % not sure what this means
BuffersPerAcquisition = 1;

%Number of trigger events you want to capture. Auto Settings do not work
%with this stage.
NumberofTriggerEvents = 1;

%End of user editable part.
%%
load('\\SimpsonNAS2\Data\SashaRazumtcev\UVF-PTIR\10-21-21 Griseofulvin spectrum attempt 2\Hadamard_matrix.mat')
folderpath1 = '\\SimpsonNAS2\Data\SashaRazumtcev\UVF-PTIR\10-21-21 Griseofulvin spectrum attempt 2\Hadamard rawdata';
folderlist1 = dir(folderpath1);
filePattern = fullfile(folderpath1, '*.mat');
theFiles = dir(filePattern);
rawimagestack = zeros(98,98,length(theFiles));
for numberofrawdata = 1:length(theFiles)
    baseFileName = ['Raw_had' num2str(numberofrawdata) '_45us_1kHz_0.94gain_GI.mat'];
    fullFileName = fullfile(theFiles(1).folder, baseFileName);
    load(fullFileName)
    number_of_pixels = (Image_V_Pixel-1)*Image_H_Pixel;
    rawdata2 = rawdata((length(rawdata)/2+1):length(rawdata));
    rawdata1 = rawdata(1:((length(rawdata)/2)-16));

    %Next, the code recontructs a fluorescence image. It is done based on a
    %reference from the piezo stage (to know precisely when the stage is
    %scanning each line). Because this alazar cards only has 2 input channels,
    %we cannot collect fluorescence intensity, QCL reference for modulation and
    %stage reference at the same time. Instead, we have precollected the stage
    %reference and found the pixels positions for the most commonly used
    %imaging conditions. The code loads 2 files for the set of imaging
    %conditions and uses it to recontruct the image. if you need to use a new
    %set of imaging conditions, run a test acquisition with stage reference in
    %channel B and calculate the "high" and "low" files using the
    %"Ramp_reconstruction" function

    RawImage = zeros(Image_V_Pixel-1,Image_H_Pixel); %preallocate an array for the image
    pixelgap = zeros(Image_V_Pixel-1,1);

    load('rampscan_300s_50deg_high_index_50kHz.mat') %Load pixel positions reference. Change the file name to match your conditions.
    load('rampscan_300s_50deg_low_index_50kHz.mat')

    %This loop divides rawdata into pixels and averages over each pixel to
    %calculate fluorescence intensity
    for i = 1:Image_V_Pixel-1
        pixelgap(i) = (high(i+1) - low(i))/Image_H_Pixel;
        for j = 1:Image_H_Pixel
            RawImage(i,j) = sum(rawdata1(1,round(low(i)+(j-1)*pixelgap):round(low(i)+(j)*pixelgap)))/length(rawdata1(1,round(low(i)+(j-1)*pixelgap):round(low(i)+(j)*pixelgap)));
        end
    end
    number_of_pixels = (Image_V_Pixel-1)*Image_H_Pixel;
    j = 1;
    k = 0;
    pixgap = round(mean(pixelgap));
    pixels_IR = zeros(number_of_pixels,pixgap+1);
    for i = 1:number_of_pixels
        pixels_IR(i,:) = rawdata1(low(j) + k*pixgap:low(j)+(k+1)*pixgap);
        k = k+1;
        if (i+1) - j*100 > 0
            j = j+1;
            k = 0;
        end
    end
    samplesperIRperiod = MidIR_period*ConfigParameters.GlobalSamplingRate;
    pixels_avg = zeros(number_of_pixels,samplesperIRperiod);

    for i=1:number_of_pixels
       for j = 1:samplesperIRperiod
           pixel_temp = 0;
           for k = 1:(pixgap/samplesperIRperiod - 1)
               pixel_temp2 = pixels_IR(i,j+(k-1)*samplesperIRperiod);
               pixel_temp = pixel_temp + pixel_temp2;
           end
           pixels_avg(i,j) = pixel_temp/(pixgap/samplesperIRperiod);
       end
    end

    %Calculate modulation depth simply by taking the difference between the
    %initial fluorescence intensity and intensity during an IR firing event

    pixels_mod = max(pixels_avg,[],2) - min(pixels_avg,[],2);

    %Generate and plot image of modulation depth at each pixel

    Image_modulation = reshape(pixels_mod,[Image_V_Pixel,Image_H_Pixel-1]);
    Image_modulation = Image_modulation';

    rawimagestack(:,:,numberofrawdata) = Image_modulation(1:98,1:98);
    
end
figure;
imagesc(rawimagestack(:,:,1));
colormap(jet);
position = getrect;
position = round(position);
SS2 = squeeze(sum(sum(rawimagestack(position(2):position(2)+position(4),position(1):position(1)+position(3),:))));
figure, plot(SS2)
% SS = squeeze(sum(sum(rawimagestack)));
% load('Z:\SashaRazumtcev\UVF-PTIR\10-19-21 Griseofulvin spectrum attempt 1\Hadamard_matrix.mat')
R = SS2'*inv(H_Mat(1:31,1:31));
figure,plot(R)

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
    baseFileName = ['rawdata_had' num2str(numberofrawdata) '.mat'];
    if numberofrawdata == 22
        baseFileName = ['rawdata_had' num2str(23) '.mat'];
    end
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