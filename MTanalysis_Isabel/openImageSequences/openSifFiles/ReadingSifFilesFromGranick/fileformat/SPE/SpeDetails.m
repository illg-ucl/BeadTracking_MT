function fileinfo=SpeDetails(file)

%fileinfo=SpeDetails(spe)
%
%SpeDetails.m is designed to scan the header of a .SPE file, as typically used
%by WinSpec, and extract the relevant parameters, including where the 
%binary data begins. After this, other programs like SpeFrame.m are then 
%able to read in frames of the movie.
%
%INCLUDE:   
%
%INPUTS:    FILE:       The path and filename of the file to be examined,
%                       ('C:\path\file.spe')
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
%   SPE Specific
%           .bitformat  The bit format of the data, whether 16 or 32 bit.
%                       There are 4 possible values:
%                           0 =   float (4 bytes)           'float32'
%                           1 =   long (4 bytes)            'int32'
%                           2 =   short (2 bytes)           'int16'
%                           3 =   unsigned short (2 bytes)  'uint16'
%           .readout:   At present, there is some uncertainty whether the
%                       time between frames is the xposure, the readout, or
%                       the sum. 
%           .pixels:    The pixels used as calibration points, to which the
%                       wavelengths were mapped 
%           .wavelengths:   The wavelengths to which the calibration pixels
%                       correspond 
%
%Copyright Stephen Anthony 2004 U. Illinois Urbana-Champaign
%Last modified by Stephen Anthony on 11/04/2009

%Open the file for reading
fid=fopen(file,'r'); 

%If clause checks if file exists.
if fid>0
    %Read in the file header, which stores information about the movie
    %dimensions. 
    header=fread(fid,2050,'uint16=>uint16');
    %;%4100bytes/2
    %value here *2 -2
    
    %Pick out the portions of the header that store the dimensions
    dimX = double(header(22));
	dimY = double(header(329));
	dimZ = double(header(724));
    
    %Determine the number of valid calibration pairs
    fseek(fid,3102,-1);
    cpairs=fread(fid,1,'uint8');
    
    %Pick out the calibration pixels (Winspec allows up to 10)
    fseek(fid,3103,-1);
    pxels=fread(fid,cpairs,'double');
        
    %Pick out the calibration wavelengths
    fseek(fid,3183,-1);
    wvlngths=fread(fid,cpairs,'double');

    %The actual camera integration time (not presently used)
%    fseek(fid,10,-1);
%    integration=fread(fid,1,'float32')*1000;

    %Pick out the time between frames in ms
    fseek(fid,672,-1);
    xposure=fread(fid,1,'float32');
    
    %Determine the format. 
    bitformat=double(header(55));
    
    fclose(fid);
else
    disp(['File ' spe ' does not exist']);
end

%Output the desired information in a generic structure format. 
fileinfo.location=file;
fileinfo.format='.spe';
fileinfo.dimX=dimX;
fileinfo.dimY=dimY;
fileinfo.Nframes=dimZ;
fileinfo.xposure=xposure;

fileinfo.readout=xposure;
fileinfo.pixels=pxels;
fileinfo.wavelengths=wvlngths;
fileinfo.bitformat=bitformat;

