function fileinfo=SifDetails(file)

%fileinfo=SifDetails(file)
%
%SifDetails.m is designed to scan the header of an Andor Technology
%Multi-Channel File (.sif file), and extract the relevant parameters,
%including where the binary data begins. After this, other programs like
%SifFrame.m are then able to read in frames of the movie.
%
%INCLUDE:   LOGICOR.m
%
%INPUTS:    FILE:       The path and filename of the file to be examined,
%                       ('C:\path\file.sif')
%
%OUTPUTS:   FILEINFO:   A structure containing pertinent information about
%                       the file. Some are universal, regardless of file
%                       format, while others are only relevant to this file
%                       format. 
%   Generic
%           .location   The complete path and file name of the file. 
%           .format     The file format, in this case 'sif'
%           .dimX,Y:    The spatial dimensions of the movie file.
%           .Nframes:   The number of frames in the movie
%           .xposure:   The total time between frames in milliseconds.
%                       (This may be longer than the exposure time.)
%   SIF Specific
%           .datastart: A position variable giving the start of data in the
%                       sif file.
%           .xyoffset:  If the absolute location of the image matters, such
%                       as for calibrated spectroscopy, this gives
%                       specifies the position on the CCD camera.
%
%Copyright Stephen Anthony 1/2009 U. Illinois Urbana-Champaign
%Last modified by Stephen Anthony on 10/03/2009


%Open up the desired file.
fid=fopen(file,'r');
if fid < 0
   error('Could not open the file.');
end

%Adjust the position to the desired positions
fseek(fid,0,0);

%Initialize the cell array to store all information from the header
variables=cell(143,2);

%Begin reading the header
variables{1,2}='File Type';
variables{1,1}=fgetl(fid);

%Confirm this is an appropriate Andor file.
if ~isequal(variables{1,1},'Andor Technology Multi-Channel File')
   fclose(fid);
   error('Not the expected Andor SIF image file.');
end

variables{2,2}='File Version';
variables(2,1)=readword(fid);
%Verify appropriate version
if ~logicor(str2double(variables{2,1}),65538)
   fclose(fid);
   error('Unexpected file version. May work, but needs verification.');
end


variables{3,2}='Signal Present';
variables(3,1)=readline(fid);

%Determine if the signal is present
if str2double(variables{3,1})==1
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TInstaImage
    variables{4,2}='version';
    variables(4,1)=readword(fid);
    %Verify appropriate version
    if ~logicor(str2double(variables{4,1}),[65550 65552 65555 65558 65559])
       fclose(fid);
       error('Unexpected InstaImage version. May work, but needs verification.');
    end

    variables{5,2}='type';
    variables(5,1)=readword(fid);
    
    variables{6,2}='active';
    variables(6,1)=readword(fid);
    
    variables{7,2}='structure_version';
    variables(7,1)=readword(fid);
    
    variables{8,2}='timedate';
    variables(8,1)=readword(fid);
    
    variables{9,2}='temperature';
    variables(9,1)=readword(fid);
    
    variables{10,2}='head';
    variables{11,2}='store_type';
    variables{12,2}='data_type';
    variables{13,2}='mode';
    variables{14,2}='trigger_source';
    skipBytes(fid,10)
    
    variables{15,2}='trigger_level';
    variables(15,1)=readword(fid);
    
    variables{16,2}='exposure_time';
    variables(16,1)=readword(fid);
    
    variables{17,2}='delay';
    variables(17,1)=readword(fid);
    
    variables{18,2}='integration_cycle_time';
    variables(18,1)=readword(fid);
    
    variables{19,2}='no_integrations';
    variables(19,1)=readword(fid);
    
    variables{20,2}='sync';
    skipBytes(fid,2)
    
    variables{21,2}='kinetic_cycle_time';
    variables(21,1)=readword(fid);
    
    variables{22,2}='pixel_readout_time';
    variables(22,1)=readword(fid);
    
    variables{23,2}='no_points';
    variables(23,1)=readword(fid);
    
    variables{24,2}='fast_track_height';
    variables(24,1)=readword(fid);
    
    variables{25,2}='gain';
    variables(25,1)=readword(fid);
    
    variables{26,2}='gate_delay';
    variables(26,1)=readword(fid);
    
    variables{27,2}='gate_width';
    variables(27,1)=readword(fid);
    
    variables{28,2}='gate_step';
    variables(28,1)=readword(fid);
    
    variables{29,2}='track_height';
    variables(29,1)=readword(fid);
    
    variables{30,2}='series_length';
    variables(30,1)=readword(fid);

    
    variables{31,2}='read_pattern';
    variables{32,2}='shutter_delay';
    skipBytes(fid,4)
    
    variables{33,2}='st_centre_row';
    variables(33,1)=readword(fid);
    
    variables{34,2}='mt_offset';
    variables(34,1)=readword(fid);
    
    variables{35,2}='operation_mode';
    variables(35,1)=readword(fid);
    
    variables{36,2}='FlipX';
    variables(36,1)=readword(fid);
    
    variables{37,2}='FlipY';
    variables(37,1)=readword(fid);
    
    variables{38,2}='Clock';
    variables(38,1)=readword(fid);
    
    variables{39,2}='AClock';
    variables(39,1)=readword(fid);
    
    variables{40,2}='MCP';
    variables(40,1)=readword(fid);
    
    variables{41,2}='Prop';
    variables(41,1)=readword(fid);
    
    variables{42,2}='IOC';
    variables(42,1)=readword(fid);
    
    variables{43,2}='Freq';
    variables(43,1)=readword(fid);
    
    variables{44,2}='VertClockAmp';
    variables(44,1)=readword(fid);
    
    variables{45,2}='data_v_shift_speed';
    variables(45,1)=readword(fid);
    
    variables{46,2}='OutputAmp';
    variables(46,1)=readword(fid);
    
    variables{47,2}='PreAmpGain';
    variables(47,1)=readword(fid);
    
    variables{48,2}='Serial';
    variables(48,1)=readword(fid);
    
    variables{49,2}='NumPulses';
    variables(49,1)=readword(fid);
    
    variables{50,2}='mFrameTransferAcqMode';
    variables(50,1)=readword(fid);
    
    if str2double(variables{4,1})~=65550
        variables{51,2}='unstabilizedTemperature';
        variables(51,1)=readword(fid);
    
        variables{52,2}='mBaselineClamp';
        variables(52,1)=readword(fid);
    end
    
    variables{53,2}='mPreScan';
    variables(53,1)=readword(fid);

    if logicor(str2double(variables{4,1}),[65555 65558])
        variables{54,2}='mEMRealGain';
        variables(54,1)=readword(fid);
    
        variables{55,2}='mBaselineOffset';
        variables(55,1)=readword(fid);
    
        variables{56,2}='mSWVersion';
        variables(56,1)=readword(fid);
    end

    if str2double(variables{4,1})==65558
        variables{57,2}='miGateMode';
        variables(57,1)=readword(fid);
    
        variables{58,2}='mSWDllVer';
        variables(58,1)=readword(fid);
    
        variables{59,2}='mSWDllRev';
        variables(59,1)=readword(fid);
    
        variables{60,2}='mSWDllRel';
        variables(60,1)=readword(fid);
    
        variables{61,2}='mSWDllBld';
        variables(61,1)=readword(fid);
    end
    
    variables{62,2}='len';
    variables(62,1)=readline(fid);

    variables{63,2}='head_model';
    variables(63,1)=readline(fid);
    
    variables{64,2}='detector_format_x';
    variables(64,1)=readword(fid);
    
    variables{65,2}='detector_format_y';
    variables(65,1)=readword(fid);
    
    variables{66,2}='len';
    variables(66,1)=readline(fid);
    
    variables{67,2}='filename';
    variables(67,1)=readline(fid);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TUserText

    variables{68,2}='TUserText Version';
    variables(68,1)=readword(fid);
    
    %Verify appropriate version
    if ~logicor(str2double(variables{68,1}),65538)
       fclose(fid);
       error('Unexpected UserText version. May work, but needs verification.');
    end

    
    variables{69,2}='len';
    variables(69,1)=readline(fid);
    
    variables{70,2}='User Text';
    variables(70,1)=readline(fid);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Tshutter
    variables{71,2}='TShutter Version';
    variables(71,1)=readword(fid);
    %Verify appropriate version
    if ~logicor(str2double(variables{68,1}),65538)
       fclose(fid);
       error('Unexpected Shutter version. May work, but needs verification.');
    end
    
    variables{72,2}='type';
    variables(72,1)=readword(fid);
    
    variables{73,2}='mode';
    variables(73,1)=readword(fid);
    
    variables{74,2}='custom_bg_mode';
    variables(74,1)=readword(fid);

    variables{75,2}='custom_mode';
    variables(75,1)=readword(fid);
    
    variables{76,2}='closing_time';
    variables(76,1)=readword(fid);

    variables{76,2}='opening_time';
    variables(76,1)=readline(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TShamrockSave
    variables{77,2}='TShamrockSave Version';
    variables(77,1)=readword(fid);
    %Verify appropriate version
    if ~logicor(str2double(variables{77,1}),[65536 65538 65540])
       fclose(fid);
       error('Unexpected ShamrockSave version. May work, but needs verification.');
    end

    variables{78,2}='isActive';
    variables(78,1)=readword(fid);

    variables{79,2}='waveDrivePresent';
    variables(79,1)=readword(fid);
    
    variables{80,2}='wavelength';
    variables(80,1)=readword(fid);
    
    variables{81,2}='gratingTurretPresent';
    variables(81,1)=readword(fid);
    
    variables{82,2}='grating';
    variables(82,1)=readword(fid);

    variables{83,2}='gratingLines';
    variables(83,1)=readword(fid);
    
    variables{84,2}='gratingBlaze';
    variables(84,1)=readline(fid);
    
    variables{85,2}='slitPresent';
    variables(85,1)=readword(fid);
    
    variables{86,2}='slitWidth';
    variables(86,1)=readword(fid);
    
    variables{87,2}='flipperMirrorPresent';
    variables(87,1)=readword(fid);
    
    variables{88,2}='flipperPort';
    variables(88,1)=readword(fid);
    
    variables{89,2}='filterPresent';
    variables(89,1)=readword(fid);
    
    variables{90,2}='filterIndex';
    variables(90,1)=readword(fid);
    
    variables{91,2}='len';
    variables(91,1)=readword(fid);
    
    variables{92,2}='filterLabel';
    variables(92,1)=readword(fid);
    
    variables{93,2}='accessoryAttached';
    variables(93,1)=readword(fid);
    
    variables{94,2}='port1State';
    variables(94,1)=readword(fid);
    
    variables{95,2}='port2State';
    variables(95,1)=readword(fid);
    
    variables{96,2}='port3State';
    variables(96,1)=readword(fid);
    
    variables{97,2}='inputPortState';
    variables(97,1)=readword(fid);
    
    variables{98,2}='outputSlitPresent';
    variables(98,1)=readword(fid);

    variables{99,2}='outputSlitWidth';
    variables(99,1)=readline(fid);


    if str2double(variables{77,1})==65540
        
        variables{100,2}='IsStepAndGlue';
        variables(100,1)=readword(fid);

        variables{101,2}='SpectrographName';
        variables(101,1)=readline(fid);
        
        %Unknown parameters
        readword(fid);
        readword(fid);
        readword(fid);
        readword(fid);
        readword(fid);
        readword(fid);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if logicor(str2double(variables{4,1}),[65558 65559])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TSpectrographSave
        variables{102,2}='TSpectrographSave Version';
        variables(102,1)=readword(fid);
        %Verify appropriate version
        if ~logicor(str2double(variables{102,1}),[65536 65539])
           fclose(fid);
           error('Unexpected TSpectrographSave version. May work, but needs verification.');
        end

        variables{103,2}='isActive';
        variables(103,1)=readword(fid);

        variables{104,2}='wavelength';
        variables(104,1)=readword(fid);

        variables{105,2}='gratingLines';
        variables(105,1)=readword(fid);

        variables{106,2}='SpectrographName';
        variables(106,1)=readline(fid);
        
        variables{107,2}='Unspecified';
        variables(107,1)=readword(fid);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TCalibImage
    variables{108,2}='TCalibImage Version';
    variables(108,1)=readword(fid);
    
    %Verify appropriate version
    if ~logicor(str2double(variables{108,1}),65539)
        %One extra line may have been added. Test this. 
        variables(108,1)=readword(fid);
        if ~logicor(str2double(variables{108,1}),65539)
            fclose(fid);
            error('Unexpected TCalibImage version. May work, but needs verification.');
        end
    end
    
    variables{109,2}='Various';
    variables(109,1)=readline(fid);
    
    variables{110,2}='x_cal0';
    variables(110,1)=readword(fid);
    
    variables{111,2}='x_cal1';
    variables(111,1)=readword(fid);

    variables{112,2}='x_cal2';
    variables(112,1)=readword(fid);
    
    variables{113,2}='x_cal3';
    variables(113,1)=readline(fid);

    variables{114,2}='y_cal0';
    variables(114,1)=readword(fid);
    
    variables{115,2}='y_cal1';
    variables(115,1)=readword(fid);

    variables{116,2}='y_cal2';
    variables(116,1)=readword(fid);
    
    variables{117,2}='y_cal3';
    variables(117,1)=readline(fid);
    
    variables{118,2}='z_cal0';
    variables(118,1)=readword(fid);
    
    variables{119,2}='z_cal1';
    variables(119,1)=readword(fid);

    variables{120,2}='z_cal2';
    variables(120,1)=readword(fid);
    
    variables{121,2}='z_cal3';
    variables(121,1)=readline(fid);
    
    variables{122,2}='rayleigh_wavelength';
    variables(122,1)=readline(fid);
    
    variables{123,2}='pixel_length';
    variables(123,1)=readline(fid);
    
    variables{124,2}='pixel_height';
    variables(124,1)=readline(fid);
    
    variables{125,2}='lenx';
    variables(125,1)=readline(fid);
    skipBytes(fid,str2double(variables{125,1}))
    
    variables{126,2}='leny';
    variables(126,1)=readline(fid);
    skipBytes(fid,str2double(variables{126,1}))

    variables{127,2}='lenz';
    variables(127,1)=readline(fid);
    skipBytes(fid,str2double(variables{127,1}))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %TImage
    variables{128,2}='TImage Version';
    variables(128,1)=readword(fid);
    %Verify appropriate version
    if ~logicor(str2double(variables{128,1}),[65538 65541])
       fclose(fid);
       error('Unexpected TImage version. May work, but needs verification.');
    end

    variables{129,2}='image_format.left';
    variables(129,1)=readword(fid);

    variables{130,2}='image_format.top';
    variables(130,1)=readword(fid);

    variables{131,2}='image_format.right';
    variables(131,1)=readword(fid);
    
    variables{132,2}='image_format.bottom';
    variables(132,1)=readword(fid);
    
    %The number of frames
    variables{133,2}='no_images';
    variables(133,1)=readword(fid);
    
    variables{134,2}='no_subimages';
    variables(134,1)=readword(fid);

    variables{135,2}='total_length';
    variables(135,1)=readword(fid);

    variables{136,2}='image_length';
    variables(136,1)=readline(fid);

    %Repeat # of subimages times
    for j=1:str2double(variables{134,1})
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TSubimage 
        variables{137,2}='TSubImage Version';
        variables(137,1)=readword(fid);
        %Verify appropriate version
        if ~logicor(str2double(variables{137,1}),65538)
           fclose(fid);
           error('Unexpected TSubImage version. May work, but needs verification.');
        end

        variables{138,2}='left';
        variables(138,1)=readword(fid);
        
        variables{139,2}='top';
        variables(139,1)=readword(fid);

        variables{140,2}='right';
        variables(140,1)=readword(fid);

        variables{141,2}='bottom';
        variables(141,1)=readword(fid);
        
        variables{142,2}='vertical_bin';
        variables(142,1)=readword(fid);

        variables{142,2}='horizontal_bin';
        variables(142,1)=readword(fid);

        variables{143,2}='subimage_offset';
        variables(143,1)=readline(fid);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %Repeat # of images times
    for j=1:str2double(variables{133,1})
        fgetl(fid);
    end
end

if logicor(str2double(variables{4,1}),65559)
   readline(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output the variables we need
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The number of frames
dimZ=str2double(variables{133,1});

%Where the binary data starts
datastart=ftell(fid);

%The number of pixels in the x and y directions
dimX=str2double(variables{140,1})-str2double(variables{138,1})+1;
dimY=str2double(variables{139,1})-str2double(variables{141,1})+1;

%The exposure time, or the time between frames, in milliseconds
xposure=str2double(variables{21,1})*1000;

%The location of the first pixel, if it matters
xyoffset=[str2double(variables{138,1}) str2double(variables{141,1})]-1;

%Close the file
fclose(fid);

%Output the desired information in a generic structure format. 
fileinfo.location=file;
fileinfo.format='.sif';
fileinfo.dimX=dimX;
fileinfo.dimY=dimY;
fileinfo.Nframes=dimZ;
fileinfo.xposure=xposure;
fileinfo.datastart=datastart;
fileinfo.xyoffset=xyoffset;



function out=readword(fid)

%Read each word, removing trailing white space
tmp=textscan(fid,'%s ',1);
out=tmp{:};

function out=readline(fid)

%Read each line, removing trailing white space
out={deblank(fgetl(fid))};


function skipBytes(fid,N)

%Skip bytes.
%
% fid      File handle
% N      Number of bytes to skip

[s,n]=fread(fid,N,'uint8');
if n < N
   fclose(fid);
   error('Inconsistent image header.');
end