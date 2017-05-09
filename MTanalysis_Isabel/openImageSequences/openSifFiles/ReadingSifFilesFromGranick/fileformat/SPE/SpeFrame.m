function arr=SpeFrame(fileinfo,frame)

%arr=SpeFrame(fileinfo,frame)
%
%This program allows frame by frame readin of the .SPE file format used by
%WinSpec. Memory limitations make it advisable to read individual frames, 
%rather than the entire movie at once, when possible. 
%
%REQUIRES:  SpeDetails
%
%INPUTS:    FILEINFO:   The output of SpeDetails.m or FileDetails.m
%           FRAME:      The frame number to read
%
%OUTPUTS:   ARR:     A matrix representation of the relevant frame.
%
%Copyright Stephen Anthony 12/2005 U. Illinois Urbana-Champaign
%Last modified by Stephen Anthony on 11/04/2009

%Vary the behavior depending upon the number of bits.
bittype=fileinfo.bitformat;

%Open up the desired file.
fid=fopen(fileinfo.location,'r'); 

%Determine the number of pixels in a frame.
Npixels=fileinfo.dimX*fileinfo.dimY;

if bittype==0 || bittype==1
    %Either 32 bit format
    bits=32;
elseif bittype==2 || bittype==3
    %Either 32 bit format
    bits=16;    
else
    error('Unsupported number of bits')
end

%Adjust the position to the beginning of the desired frame. 
newpos=4100+Npixels*(frame-1)*bits/8;
fseek(fid,newpos,0);

if bittype==0
    %Read in the 32-bit float frame
    arr=double(fread(fid,Npixels,'float32=>float32'));
elseif bittype==1
    %Read in the 32-bit signed integer frame
    arr=double(fread(fid,Npixels,'int32=>int32'));    
elseif bittype==2
    %Read in the 16-bit signed integer frame
    arr=double(fread(fid,Npixels,'int16=>int16'));
elseif bittype==3
    %Read in the 16-bit unsigned integer frame
    arr=double(fread(fid,Npixels,'uint16=>uint16'));
end

%Convert it into a matrix, and then transpose it to orient it normally.
arr=reshape(arr,fileinfo.dimX,fileinfo.dimY)';

%Close the file
fclose(fid);

%The following information is not needed for this program, but is included to
%as helpful for anyone needing to debug this. 
%
%position=ftell(fid)
%
%This command tells the current position of reading in. The file format for
%spe contains a 4100 byte header, followed by the actual data. Assuming a
%movie of N frames with dimensions X*Y, the first X values are the first
%row, the next X values are the second row, etc. The first X*Y values
%therefore correspond to the first frame, the second X*Y values to the
%second frame, etc. Note, each value is stored as a 16-bit integer, or two
%bytes, so the position of the 27th frame would be 4100+26*X*Y*2, the 4100
%corresponding to the offset for the header. 
%
%Fopen opens the file at the beginning. after running fread, the position
%is now wherever you stopped reading. Ftell tells where you currently are.
%Fseek is used to manually change your position. This program uses Fseek to
%position the reading start point at the beginning of the desired frame,
%and then read X*Y numbers in, corresponding to the frame, followed by
%converting them from a column to a matrix of the orientation seen in
%Winspec. 
