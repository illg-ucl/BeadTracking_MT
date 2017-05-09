function fileinfo=FileDetails(file)

%fileinfo=FileDetails(file)
%
%FileDetails.m is a front end which unifies the commands necessary to scan
%file headers and determine relevant parameters. 
%
%INCLUDE:   LOGICOR.m
%           SpeDetails.m
%           SifDetails.m
%           CineDetails.m
%           RootDetails.m
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
%           .format     The file format, in this case '.ext'
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
%   SPE Specific
%           .bitformat  Whether this is a 16-bit or 32-bit SPE
%           .readout:   At present, there is some uncertainty whether the
%                       time between frames is the xposure, the readout, or
%                       the sum. 
%           .pixels:    The pixels used as calibration points, to which the
%                       wavelengths were mapped 
%           .wavelengths:   The wavelengths to which the calibration pixels
%                       correspond
%
%Copyright Stephen Anthony 10/2009 U. Illinois Urbana-Champaign
%Last modified by Stephen Anthony on 11/04/2009

%Determine the type of file based upon extension
[directory,name,ext]=fileparts(file);

%Call the appropriate function dependent upon file type, by comparing the
%file extension in a non-case sensitive fashion
if strcmpi(ext,'.sif')
    %Andor Technology MultiChannel File Format
    fileinfo=SifDetails(file);
elseif strcmpi(ext,'.spe')
    %SPE File format (WinSpec)
    fileinfo=SpeDetails(file);
elseif strcmpi(ext,'.cine')
    %CINE Vision Research File Format
    fileinfo=CineDetails(file);
elseif strcmpi(ext,'.root')
    %ROOT File Format
    fileinfo=RootDetails(file);
end