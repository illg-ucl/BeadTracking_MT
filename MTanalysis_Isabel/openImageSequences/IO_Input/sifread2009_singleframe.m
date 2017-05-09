function [newdata,headStr] = sifread2009_singleframe(file,whichframe,bright_field_image)
% 
% read SINGLE frame from *.sif file
%
% file -- file name
% whichframe -- frame number (index begins from 1 instead of 0)
% bright_field_image -- flag to point out the current image is bright field image 
%
% there is a bug 
%   bright field image should be rotated by 90 degree to consist with the
%   result from Andor softare
%
% newdata -- current frame image data
% headStr -- head information
%

% % [ReturnCode, SeriesLength, ImageSize, TotalAcquisitionSize] = GetAndorSifSize(file,0);
% % % fixing bug!
% % % if whichframe>=SeriesLength
% % %     return;
% % % end
% % 
% % 
% % % Set file access mode to read all
% % access_mode = 0;
% % rc=AndorSifSetAccessMode(access_mode);
% % 
% % % read the file from specific location ** cHANGE THIS **
% % 
% % rc=AndorSifReadFromFile(file);
% % 
% % % check it is loaded
% % [rc,loaded]=AndorSifIsLoaded;
% % 
% % signal=0;
% % % check that signal is present
% % [rc,present]=AndorSifIsDataSourcePresent(signal);
% % 
% % % Get the number of signal frames in the sif
% % [rc,no_frames]=AndorSifGetNumberFrames(signal);
% % 
% % % Determine how many sub images make up 1 frame
% % [rc,no_subImgs]=AndorSifGetNumberSubImages(signal);
% % 
% % % get the sub image info (zero based index);
% % [rc,left,bottom,right,top,vBin,hBin]=AndorSifGetSubImageInfo(signal,0);
% % 
% % %  Determine the size of 1 frame
% % [rc,frameSize]=AndorSifGetFrameSize(signal);
% % 
% % % Retrieve a frame (zero based index);
% % N = whichframe; % Nth frame in the sequenc and ZERO is the first one from this sequence
% % [rc,data]=AndorSifGetFrame(signal,N,frameSize);
% % 
% % % remap the data
% % width=(right-left)+1;
% % % disp(width);
% % height=(top-bottom)+1;
% % % disp(height);
% % 
% % newdata=reshape(data,width,height);
% % if isequal(bright_field_image,1)
% %     newdata = imrotate(newdata,90); % Must done because of consistent with Andor (quan is confusing about this!)
% % end
% %     
% % 
% % % imwidth  = sqrt(ImageSize);
% % % imheight = sqrt(ImageSize);
% % % pointerStart  = 1;
% % % pointerEnd = 0;
% % % for iframe=1:round(SeriesLength) % this line is only used for debug
% % % % for iframe=1:SeriesLength
% % %     pointerStart = 1+(iframe-1)*imwidth*imheight;
% % %     pointerEnd = pointerStart-1+imwidth*imheight;
% % %     oneFrame1D = ResultArray(pointerStart:pointerEnd);
% % %     oneFrame2D = reshape(oneFrame1D,[imwidth imheight]);
% % %     oneFrame = imrotate(oneFrame2D,90,'bicubic');
% % %     dataArry{iframe}.image = oneFrame;
% % % end
% % 
% %     [ReturnCode, headStr.AClock] = GetAndorSifProperty(file, 'AClock', 0);
% %     [ReturnCode, headStr.Active] = GetAndorSifProperty(file, 'Active', 0);
% %     [ReturnCode, headStr.BaselineClamp] = GetAndorSifProperty(file, 'BaselineClamp', 0);
% %     [ReturnCode, headStr.BaselineOffset] = GetAndorSifProperty(file, 'BaselineOffset', 0);
% %     [ReturnCode, headStr.CentreRow] = GetAndorSifProperty(file, 'CentreRow', 0);
% %     [ReturnCode, headStr.Clock] = GetAndorSifProperty(file, 'Clock', 0);
% %     [ReturnCode, headStr.DataType] = GetAndorSifProperty(file, 'DataType', 0);
% %     [ReturnCode, headStr.Delay] = GetAndorSifProperty(file, 'Delay', 0);
% %     [ReturnCode, headStr.DetectorFormatX] = GetAndorSifProperty(file, 'DetectorFormatX', 0);
% %     [ReturnCode, headStr.DetectorFormatZ] = GetAndorSifProperty(file, 'DetectorFormatZ', 0);
% %     [ReturnCode, headStr.EMRealGain] = GetAndorSifProperty(file, 'EMRealGain', 0);
% %     [ReturnCode, headStr.ExposureTime] = GetAndorSifProperty(file, 'ExposureTime', 0);
% %     [ReturnCode, headStr.FileName] = GetAndorSifProperty(file, 'FileName', 0);
% %     [ReturnCode, headStr.FlipX] = GetAndorSifProperty(file, 'FlipX', 0);
% %     [ReturnCode, headStr.FlipY] = GetAndorSifProperty(file, 'FlipY', 0);
% %     [ReturnCode, headStr.FormattedTime] = GetAndorSifProperty(file, 'FormattedTime', 0);
% %     [ReturnCode, headStr.FrameTransferAcquisitionMode] = GetAndorSifProperty(file, 'FrameTransferAcquisitionMode', 0);
% %     [ReturnCode, headStr.Frequency] = GetAndorSifProperty(file, 'Frequency', 0);
% %     [ReturnCode, headStr.Gain] = GetAndorSifProperty(file, 'Gain', 0);
% %     [ReturnCode, headStr.Head] = GetAndorSifProperty(file, 'Head', 0);
% %     [ReturnCode, headStr.HeadModel] = GetAndorSifProperty(file, 'HeadModel', 0);
% %     [ReturnCode, headStr.IntegrationCycleTime] = GetAndorSifProperty(file, 'IntegrationCycleTime', 0);
% %     [ReturnCode, headStr.IOC] = GetAndorSifProperty(file, 'IOC', 0);
% %     [ReturnCode, headStr.KineticCycleTime] = GetAndorSifProperty(file, 'KineticCycleTime', 0);
% %     [ReturnCode, headStr.MCP] = GetAndorSifProperty(file, 'MCP', 0);
% %     [ReturnCode, headStr.Mode] = GetAndorSifProperty(file, 'Mode', 0);
% %     [ReturnCode, headStr.ModeFullName] = GetAndorSifProperty(file, 'ModeFullName', 0);
% %     [ReturnCode, headStr.NumberImages] = GetAndorSifProperty(file, 'NumberImages', 0);
% %     [ReturnCode, headStr.NumberIntegrations] = GetAndorSifProperty(file, 'NumberIntegrations', 0);
% %     [ReturnCode, headStr.NumberPulses] = GetAndorSifProperty(file, 'NumberPulses', 0);
% %     [ReturnCode, headStr.Operation] = GetAndorSifProperty(file, 'Operation', 0);
% %     [ReturnCode, headStr.OutputAmplifier] = GetAndorSifProperty(file, 'OutputAmplifier', 0);
% %     [ReturnCode, headStr.PixelReadOutTime] = GetAndorSifProperty(file, 'PixelReadOutTime', 0);
% %     [ReturnCode, headStr.PreAmplifierGain] = GetAndorSifProperty(file, 'PreAmplifierGain', 0);
% %     [ReturnCode, headStr.PreScan] = GetAndorSifProperty(file, 'PreScan', 0);
% %     [ReturnCode, headStr.ReadPattern] = GetAndorSifProperty(file, 'ReadPattern', 0);
% %     [ReturnCode, headStr.ReadPatternFullName] = GetAndorSifProperty(file, 'ReadPatternFullName', 0);
% %     [ReturnCode, headStr.RowOffset] = GetAndorSifProperty(file, 'RowOffset', 0);
% %     [ReturnCode, headStr.Serial] = GetAndorSifProperty(file, 'Serial', 0);
% %     [ReturnCode, headStr.ShutterDelay] = GetAndorSifProperty(file, 'ShutterDelay', 0);
% %     [ReturnCode, headStr.SIDisplacement] = GetAndorSifProperty(file, 'SIDisplacement', 0);
% %     [ReturnCode, headStr.SINumberSubFrames] = GetAndorSifProperty(file, 'SINumberSubFrames', 0);
% %     [ReturnCode, headStr.StoreType] = GetAndorSifProperty(file, 'StoreType', 0);
% %     [ReturnCode, headStr.SWVersion] = GetAndorSifProperty(file, 'SWVersion', 0);
% %     [ReturnCode, headStr.SWVersionEx] = GetAndorSifProperty(file, 'SWVersionEx', 0);
% %     [ReturnCode, headStr.Temperature] = GetAndorSifProperty(file, 'Temperature', 0);
% %     [ReturnCode, headStr.Time] = GetAndorSifProperty(file, 'Time', 0);
% %     [ReturnCode, headStr.TrackHeight] = GetAndorSifProperty(file, 'TrackHeight', 0);
% %     [ReturnCode, headStr.TriggerLevel] = GetAndorSifProperty(file, 'TriggerLevel', 0);
% %     [ReturnCode, headStr.TriggerSource] = GetAndorSifProperty(file, 'TriggerSource', 0);
% %     [ReturnCode, headStr.TriggerSourceFullName] = GetAndorSifProperty(file, 'TriggerSourceFullName', 0);
% %     [ReturnCode, headStr.Type] = GetAndorSifProperty(file, 'Type', 0);
% %     [ReturnCode, headStr.UnstabalizedTemperature] = GetAndorSifProperty(file, 'UnstabalizedTemperature', 0);
% %     [ReturnCode, headStr.Version] = GetAndorSifProperty(file, 'Version', 0);
% %     [ReturnCode, headStr.VerticalClockAmp] = GetAndorSifProperty(file, 'VerticalClockAmp', 0);
% %     [ReturnCode, headStr.VerticalShiftSpeed] = GetAndorSifProperty(file, 'VerticalShiftSpeed', 0);
% %     
   
    % example from Anjana
file_name_sif = file;

rc=atsif_setfileaccessmode(0);
rc=atsif_readfromfile(file_name_sif);

if (rc == 22002)
  signal=0;
  [rc,present]=atsif_isdatasourcepresent(signal);
  if present
    [rc,no_frames]=atsif_getnumberframes(signal);
    if (no_frames > 0)
        [rc,size]=atsif_getframesize(signal);
        [rc,left,bottom,right,top,hBin,vBin]=atsif_getsubimageinfo(signal,0);
        xaxis=0;
        [rc,pattern]=atsif_getpropertyvalue(signal,'ReadPattern');
        
        if(pattern == '0')
            
           calibvals = zeros(1,size);
           for i=1:size,[rc,calibvals(i)]=atsif_getpixelcalibration(signal,xaxis,(i)); 
           end
           [rc,data]=atsif_getframe(signal,0,size);
           plot(calibvals,data);      
           title('spectrum');
           [rc,xtype]=atsif_getpropertyvalue(signal,'XAxisType');
           [rc,xunit]=atsif_getpropertyvalue(signal,'XAxisUnit');
           [rc,ytype]=atsif_getpropertyvalue(signal,'YAxisType');
           [rc,yunit]=atsif_getpropertyvalue(signal,'YAxisUnit');
           xlabel({xtype;xunit});
           ylabel({ytype;yunit});
           
        elseif(pattern == '4')
            
            calibvals = zeros(1,size);
            for i=1:size
                [rc,calibvals(i)]=atsif_getpixelcalibration(signal,xaxis,(i));
            end
%             figure(1),
%             for iframe = 0:(no_frames-1)
%                 [rc,data]=atsif_getframe(signal,iframe,size);
%                 B = reshape(data,[sqrt(size) sqrt(size)]);
%                 imshow(B,[])
%                 disp(strcat('current frame is:  ',num2str(iframe)))
%                 pause(0.3)                
%             end
%             close(figure(1));
            iframe = whichframe;
            [rc,data]=atsif_getframe(signal,iframe,size);
            B = reshape(data,[sqrt(size) sqrt(size)]);
            newdata = B;

 
%             title('image');
            [rc,xtype]=atsif_getpropertyvalue(signal,'XAxisType');
            [rc,xunit]=atsif_getpropertyvalue(signal,'XAxisUnit');
            [rc,ytype]=atsif_getpropertyvalue(signal,'YAxisType');
            [rc,yunit]=atsif_getpropertyvalue(signal,'YAxisUnit');
            [rc,ExposureTime]=atsif_getpropertyvalue(signal,'ExposureTime');
%             xlabel({xtype;xunit});
%             ylabel({ytype;yunit});
            headStr.DetectorFormatX = right - left   + 1;
            headStr.DetectorFormatZ = top   - bottom + 1;
            headStr.NumberImages = no_frames;
            headStr.ExposureTime = str2num(ExposureTime); % 0.003!?
        end
        
    end    
  end
  atsif_closefile;
else
  disp('Could not load file.  ERROR - ');
  disp(rc);
end
    

end % this function