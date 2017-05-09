function arr=SifFrame(fileinfo,frame)

%arr=SifFrame(fileinfo,frame)
%
%This program allows frame by frame readin of Andor's .sif file format.
%Memory limitations make it advisable to read individual frames, rather
%than the entire movie at once, when possible. 
%
%REQUIRES:  SifDetails
%
%INPUTS:    FILEINFO:   The output of SifDetails.m or FileDetails.m
%           FRAME:      The frame number to read
%
%OUTPUTS:   ARR:     A matrix representation of the relevant frame.
%
%Copyright Stephen Anthony 12/2005 U. Illinois Urbana-Champaign
%Last modified by Stephen Anthony on 10/03/2009

%Open up the desired file (with reading permission).
fid=fopen(fileinfo.location,'r'); 

%Determine the number of pixels in a frame.
Npixels=fileinfo.dimX*fileinfo.dimY;

%Adjust the position to the beginning of the desired frame. 
newpos=fileinfo.datastart+Npixels*(frame-1)*4;

%Adjust the position to the desired positions
fseek(fid,newpos,0);

%Read in the frame
arr=fread(fid,Npixels,'float32=>float32');

%Convert it into a matrix, and then transpose it to orient it normally. 
arr=reshape(arr,fileinfo.dimX,fileinfo.dimY)';
arr=double(arr);

%Close the file
fclose(fid);