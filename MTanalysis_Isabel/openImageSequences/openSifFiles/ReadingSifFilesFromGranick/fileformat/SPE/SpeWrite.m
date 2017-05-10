function SpeWrite(fileinfo,framenum,newframe)

%SifWrite(fileinfo,framenum,newframe)
%
%SifWrite.m overwrites an existing frame in a AndorFile.  This should only
%   be performed on files which have first been coppied.
%
%INCLUDE:   (none)
%
%INPUTS:    FILEINFO:   A structure containing pertinent information about
%                       the file, generated by FileDetails.m  In order for
%                       this code to work, FILEINFO must have a field .copy
%           FRAMENUM:   The index of the frame to replace
%           NEWFRAME:   The data for a new frame, where (1,1) corisponds to
%                       a pixel in the upper left corner
%
%OUTPUTS:   (none)
%               
%Copyright Scott Parker 10/2009 U. Illinois Urbana-Champaign
%Last modified by Scott Parker on 10/08/2009


% Check if the copy field exists.  If not, display a warning and terminate
% the program
if ~strcmp(fileinfo.copy,'Duplicate File')
  error('This program should only be run on copies, never on the original data')
end

bits=fileinfo.bits;

%Open up the desired file.
fid=fopen(fileinfo.location,'r+'); 

%Determine the number of pixels in a frame.
Npixels=fileinfo.dimX*fileinfo.dimY;

%Adjust the position to the beginning of the desired frame. 
newpos=4100+Npixels*(framenum-1)*bits/8;

%Adjust the position to the desired positions
fseek(fid,newpos,0);

%Adjust the image to be written so that it is all positive, and uint16.
frameZ=newframe-median(median(newframe))+1000;
%row(row>(2^31))=row(row>(2^31))-2^32;

%Reshape the frame from the image array format to the file format.
frameZ=frameZ';
frameZ=reshape(frameZ,fileinfo.dimX*fileinfo.dimY,1);


if bits==16
    frameZ=uint16(frameZ);
    fwrite(fid,frameZ, 'integer*2');
elseif bits==32
    frameZ=uint32(frameZ);
    fwrite(fid,frameZ,'integer*4');
else
    error('Unsupported number of bits')
end

%Close the file
fclose(fid);