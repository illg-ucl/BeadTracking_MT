% Read Andor SIF image files for the new version 
% Solis version: 4.12.30002.0
% SDK version: 2.84.30002.0
% Copyright@2008
%
% Synopsis:
% 
% [dataArry,headStr] = sifread2009(file)
% read sif file by name for characteristic. return image
% data and its head information in named structured as follows:
% 
%     AClock
%     Active
%     BaselineClamp
%     BaselineOffset
%     CentreRow
%     Clock
%     DataType
%     Delay
%     DetectorFormatX
%     DetectorFormatZ
%     EMRealGain
%     ExposureTime
%     FileName
%     FlipX
%     FlipY
%     FormattedTime
%     FrameTransferAcquisitionMode
%     Frequency
%     Gain
%     Head
%     HeadModel
%     IntegrationCycleTime
%     IOC
%     KineticCycleTime
%     MCP
%     Mode
%     ModeFullName
%     NumberImages
%     NumberIntegrations
%     NumberPulses
%     Operation
%     OutputAmplifier
%     PixelReadOutTime
%     PreAmplifierGain
%     PreScan
%     ReadPattern
%     ReadPatternFullName
%     RowOffset
%     Serial
%     ShutterDelay
%     SIDisplacement
%     SINumberSubFrames
%     StoreType
%     SWVersion
%     SWVersionEx
%     Temperature
%     Time
%     TrackHeight
%     TriggerLevel
%     TriggerSource
%     TriggerSourceFullName
%     Type
%     UnstabalizedTemperature
%     Version
%     VerticalClockAmp
%     VerticalShiftSpeed
% 
% Note:
%
% this function is supported by three dll provided by Andor company
%
% Quan xue @ Clarendon 14/082009
% 

function [dataArry,headStr] = sifread2009(file,bShowHead)

if nargin<1
    disp('please select one file name as inputting parameter')
    return;
elseif nargin<2
    bShowHead = 0;    
elseif nargin<3
    seq_factor = 1;
end

bReadAll = 1;

[ReturnCode, SeriesLength, ImageSize, TotalAcquisitionSize] = GetAndorSifSize(file,0);
if bReadAll == 1
    
    % ResultArray = zeros(1, TotalAcquisitionSize);
    % [ReturnCode, ResultArray] = GetAndorSifData(TotalAcquisitionSize, 0, file);
    [ReturnCode, ResultArray] = GetAndorSifData(SeriesLength*ImageSize, 0, file);
    
else
% Set file access mode to read all
access_mode = 0;
rc=AndorSifSetAccessMode(access_mode);

% read the file from specific location ** cHANGE THIS **

rc=AndorSifReadFromFile(file);

% check it is loaded
[rc,loaded]=AndorSifIsLoaded;

signal=0;
% check that signal is present
[rc,present]=AndorSifIsDataSourcePresent(signal);

% Get the number of signal frames in the sif
[rc,no_frames]=AndorSifGetNumberFrames(signal);

% Determine how many sub images make up 1 frame
[rc,no_subImgs]=AndorSifGetNumberSubImages(signal);

% get the sub image info (zero based index);
[rc,left,bottom,right,top,vBin,hBin]=AndorSifGetSubImageInfo(signal,0);

%  Determine the size of 1 frame
[rc,frameSize]=AndorSifGetFrameSize(signal);

% Retrieve a frame (zero based index);
N = 0; % Nth frame in the sequenc and ZERO is the first one from this sequence
[rc,data]=AndorSifGetFrame(signal,N,frameSize);

% remap the data
width=(right-left)+1;
% disp(width);
height=(top-bottom)+1;
% disp(height);

colormap(gray);
newdata=reshape(data,width,height);
imshow(newdata,[]);
    
end

imwidth  = sqrt(ImageSize);
imheight = sqrt(ImageSize);
pointerStart  = 1;
pointerEnd = 0;
for iframe=1:round(SeriesLength) % this line is only used for debug
% for iframe=1:SeriesLength
    pointerStart = 1+(iframe-1)*imwidth*imheight;
    pointerEnd = pointerStart-1+imwidth*imheight;
    oneFrame1D = ResultArray(pointerStart:pointerEnd);
    oneFrame2D = reshape(oneFrame1D,[imwidth imheight]);
%     oneFrame = imrotate(oneFrame2D,90,'bicubic');
    oneFrame = oneFrame2D;
    dataArry{iframe}.image = oneFrame;
end

% if bShowHead==1
    [ReturnCode, headStr.AClock] = GetAndorSifProperty(file, 'AClock', 0);
    [ReturnCode, headStr.Active] = GetAndorSifProperty(file, 'Active', 0);
    [ReturnCode, headStr.BaselineClamp] = GetAndorSifProperty(file, 'BaselineClamp', 0);
    [ReturnCode, headStr.BaselineOffset] = GetAndorSifProperty(file, 'BaselineOffset', 0);
    [ReturnCode, headStr.CentreRow] = GetAndorSifProperty(file, 'CentreRow', 0);
    [ReturnCode, headStr.Clock] = GetAndorSifProperty(file, 'Clock', 0);
    [ReturnCode, headStr.DataType] = GetAndorSifProperty(file, 'DataType', 0);
    [ReturnCode, headStr.Delay] = GetAndorSifProperty(file, 'Delay', 0);
    [ReturnCode, headStr.DetectorFormatX] = GetAndorSifProperty(file, 'DetectorFormatX', 0);
    [ReturnCode, headStr.DetectorFormatZ] = GetAndorSifProperty(file, 'DetectorFormatZ', 0);
    [ReturnCode, headStr.EMRealGain] = GetAndorSifProperty(file, 'EMRealGain', 0);
    [ReturnCode, headStr.ExposureTime] = GetAndorSifProperty(file, 'ExposureTime', 0);
    [ReturnCode, headStr.FileName] = GetAndorSifProperty(file, 'FileName', 0);
    [ReturnCode, headStr.FlipX] = GetAndorSifProperty(file, 'FlipX', 0);
    [ReturnCode, headStr.FlipY] = GetAndorSifProperty(file, 'FlipY', 0);
    [ReturnCode, headStr.FormattedTime] = GetAndorSifProperty(file, 'FormattedTime', 0);
    [ReturnCode, headStr.FrameTransferAcquisitionMode] = GetAndorSifProperty(file, 'FrameTransferAcquisitionMode', 0);
    [ReturnCode, headStr.Frequency] = GetAndorSifProperty(file, 'Frequency', 0);
    [ReturnCode, headStr.Gain] = GetAndorSifProperty(file, 'Gain', 0);
    [ReturnCode, headStr.Head] = GetAndorSifProperty(file, 'Head', 0);
    [ReturnCode, headStr.HeadModel] = GetAndorSifProperty(file, 'HeadModel', 0);
    [ReturnCode, headStr.IntegrationCycleTime] = GetAndorSifProperty(file, 'IntegrationCycleTime', 0);
    [ReturnCode, headStr.IOC] = GetAndorSifProperty(file, 'IOC', 0);
    [ReturnCode, headStr.KineticCycleTime] = GetAndorSifProperty(file, 'KineticCycleTime', 0);
    [ReturnCode, headStr.MCP] = GetAndorSifProperty(file, 'MCP', 0);
    [ReturnCode, headStr.Mode] = GetAndorSifProperty(file, 'Mode', 0);
    [ReturnCode, headStr.ModeFullName] = GetAndorSifProperty(file, 'ModeFullName', 0);
    [ReturnCode, headStr.NumberImages] = GetAndorSifProperty(file, 'NumberImages', 0);
    [ReturnCode, headStr.NumberIntegrations] = GetAndorSifProperty(file, 'NumberIntegrations', 0);
    [ReturnCode, headStr.NumberPulses] = GetAndorSifProperty(file, 'NumberPulses', 0);
    [ReturnCode, headStr.Operation] = GetAndorSifProperty(file, 'Operation', 0);
    [ReturnCode, headStr.OutputAmplifier] = GetAndorSifProperty(file, 'OutputAmplifier', 0);
    [ReturnCode, headStr.PixelReadOutTime] = GetAndorSifProperty(file, 'PixelReadOutTime', 0);
    [ReturnCode, headStr.PreAmplifierGain] = GetAndorSifProperty(file, 'PreAmplifierGain', 0);
    [ReturnCode, headStr.PreScan] = GetAndorSifProperty(file, 'PreScan', 0);
    [ReturnCode, headStr.ReadPattern] = GetAndorSifProperty(file, 'ReadPattern', 0);
    [ReturnCode, headStr.ReadPatternFullName] = GetAndorSifProperty(file, 'ReadPatternFullName', 0);
    [ReturnCode, headStr.RowOffset] = GetAndorSifProperty(file, 'RowOffset', 0);
    [ReturnCode, headStr.Serial] = GetAndorSifProperty(file, 'Serial', 0);
    [ReturnCode, headStr.ShutterDelay] = GetAndorSifProperty(file, 'ShutterDelay', 0);
    [ReturnCode, headStr.SIDisplacement] = GetAndorSifProperty(file, 'SIDisplacement', 0);
    [ReturnCode, headStr.SINumberSubFrames] = GetAndorSifProperty(file, 'SINumberSubFrames', 0);
    [ReturnCode, headStr.StoreType] = GetAndorSifProperty(file, 'StoreType', 0);
    [ReturnCode, headStr.SWVersion] = GetAndorSifProperty(file, 'SWVersion', 0);
    [ReturnCode, headStr.SWVersionEx] = GetAndorSifProperty(file, 'SWVersionEx', 0);
    [ReturnCode, headStr.Temperature] = GetAndorSifProperty(file, 'Temperature', 0);
    [ReturnCode, headStr.Time] = GetAndorSifProperty(file, 'Time', 0);
    [ReturnCode, headStr.TrackHeight] = GetAndorSifProperty(file, 'TrackHeight', 0);
    [ReturnCode, headStr.TriggerLevel] = GetAndorSifProperty(file, 'TriggerLevel', 0);
    [ReturnCode, headStr.TriggerSource] = GetAndorSifProperty(file, 'TriggerSource', 0);
    [ReturnCode, headStr.TriggerSourceFullName] = GetAndorSifProperty(file, 'TriggerSourceFullName', 0);
    [ReturnCode, headStr.Type] = GetAndorSifProperty(file, 'Type', 0);
    [ReturnCode, headStr.UnstabalizedTemperature] = GetAndorSifProperty(file, 'UnstabalizedTemperature', 0);
    [ReturnCode, headStr.Version] = GetAndorSifProperty(file, 'Version', 0);
    [ReturnCode, headStr.VerticalClockAmp] = GetAndorSifProperty(file, 'VerticalClockAmp', 0);
    [ReturnCode, headStr.VerticalShiftSpeed] = GetAndorSifProperty(file, 'VerticalShiftSpeed', 0);
% else
%     headStr = 0;
end