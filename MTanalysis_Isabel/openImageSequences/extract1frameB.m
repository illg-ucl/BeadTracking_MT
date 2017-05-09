function [frame]=extract1frameB(frame_number)
%
% ========================================
% BeadTracking_MT.
% Copyright (c) 2017. Isabel Llorente-Garcia, Dept. of Physics and Astronomy, University College London, United Kingdom.
% Released and licensed under a BSD 2-Clause License:
% https://github.com/illg-ucl/BeadTracking_MT/blob/master/LICENSE
% This program is free software: you can redistribute it and/or modify it under the terms of the BSD 2-Clause License.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the BSD 2-Clause License for more details. You should have received a copy of the BSD 2-Clause License along with this program.
% Citation: If you use this software for your data analysis please acknowledge it in your publications and cite as follows.
% -	Citation example 1: 
% BeadTracking_MT software. (Version). 2017. Isabel Llorente-Garcia, 
% Dept. of Physics and Astronomy, University College London, United Kingdom.
% https://github.com/illg-ucl/BeadTracking_MT. (Download date).
% 
% -	Citation example 2:
% @Manual{... ,
% title  = {BeadTracking_MT software. (Version).},
% author       = {{Isabel Llorente-Garcia}},
% organization = { Dept. of Physics and Astronomy, University College London, United Kingdom.},
% address      = {Gower Place, London, UK.},
% year         = 2017,
% url          = {https://github.com/illg-ucl/BeadTracking_MT}}
% ========================================
%
% This function chooses an image sequence (.sif, .dv, .tif or .mat data) 
% from an input dialog box, returns the image data corresponding to 
% "frame_number", and plots it. Based on extract_image_sequence_data.m.
% Note that it reads the full images (which is unnecessary perhaps).
%
% Example of how to call this function: extract1frame(100). 

disp('Choose input image data sequence (.sif, .dv, .tif or .mat data):')
disp(' ') % empty line

% uigetfile opens a file dialog box to choose image data file:
%[file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
[file_data,path_data] = uigetfile({'*.sif';'*.dv';'*.mat';'*.tif'},'Chose image data sequence (.sif, .dv, .tif or .mat data):');

% Error control if neither a .sif, .tif, .mat or .dv file have been selected:
if isempty(path_data) % If there is no .sif, .tif or .dv image sequence file for such image_label, show error and exit function:
    error('Check you are in the correct directory and run again. No .sif, .tif, .dv or .mat file selected.');
end

data_folder_path = strcat(path_data,file_data);
disp(data_folder_path) % write image path to command window.

data_struct = dir(data_folder_path);
% This produces something like:        
%         name: '210217r25.tif'
%        date: '29-Mar-2017 09:42:44'
%       bytes: 306759224
%       isdir: 0
%     datenum: 7.3678e+05
image_path = data_struct.name;

dvImageExists=strfind(image_path,'.dv');
sifImageExists=strfind(image_path,'.sif');
tifImageExists=strfind(image_path,'.tif');
matImageExists=strfind(image_path,'.mat');

% Turn off image size adjust warning: "Warning: Image is too big to fit on screen; displaying
% at X%". Warning identifier  is 'images:initSize:adjustingMag'.
warning('off','images:initSize:adjustingMag')


% Read in entire image first:

%% For .dv files 
% This needs bfopen.m and loci_tools.jar (see folder openDVfiles):

if isempty(dvImageExists)==0

    data = bfopen(image_path); % data is a cell array of cell arrays.
    % Assumming there is only one series, i.e., that numSeries = size(data,1)=1, we do:
    series1 = data{1,1}; % this cell array contains an array where each row (first index, p) corresponds to a frame.
    % First column (series1{p,1}) is the matrix data for the image frame. Second column (series1{p,2}) is the
    % frame label with path and info.
    series1_numFrames = size(series1,1); % number of frames in the image sequence.
    % Each image frame matrix data is series1{p,1}, with real intensity numbers
    % (not between 0 and 1).
     
    % To produce image data in final output form:
    % Loop through frames:
    for p = 1:series1_numFrames
        frame = series1{p,1};
        % frame = double(frame);  % to class double.
        image_data(p).frame_data = frame; 
    end
    
    % Add other useful info to final output:
    numFrames = series1_numFrames;
    frame_Ysize = size(image_data(1).frame_data,1);
    frame_Xsize = size(image_data(1).frame_data,2);

end


%% For .sif files:

if isempty(sifImageExists)==0
  
    % First get size of .sif image file: (see IO_Input folder, SifFunctions.txt, or page 95 my notebook 1).
    % See Isabel/myMatlabFiles/IO_input folder.
    [ReturnCode, numFrames, ImageSize, TotalAcquisitionSize]=GetAndorSifSize(image_path,0);
    
    % numFrames is the length of the sequence, the number of frames.
    % ImageSize is the size of the image, e.g. 512*512 = 262144.
    % TotalAcquisitionSize is numFrames*ImageSize
    
    % read .sif image data. sifData is an array of 1x1 structures, as many
    % columns as frames on image sequence. Read frames 1 to numFrames:
    [sifData] = read_sif_data_direct(image_path,numFrames,ImageSize,1,numFrames); % sifData is a cell array.
    
    for p = 1:numFrames
        frame = sifData{p}.sliceData; % extract frame data which is stored in the field 'sliceData'.
        % frame = double(frame); % to class double.
        image_data(p).frame_data = imrotate(frame,90); % original intensity values (needs to be rotated).
    end
    
    % Get frame dimensions:
    frame_Ysize = size(image_data(1).frame_data,1);
    frame_Xsize = size(image_data(1).frame_data,2);
    
end


%% For .tif files 

if isempty(tifImageExists)==0
    
    tif_info = imfinfo(image_path);
    % Add other useful info to final output:
    numFrames = length(tif_info); % number of frames in sequence.
    frame_Ysize = tif_info(1).Height;
    frame_Xsize = tif_info(1).Width;
%     % Make the image be a square image:
%     frame_Ysize = min(frame_Ysize,frame_Xsize);
%     frame_Xsize = min(frame_Ysize,frame_Xsize);

    % If by any chance the frame size is an odd number:
    if mod(frame_Ysize,2)~=0 % modulus after division
        frame_Ysize = frame_Ysize - 1;
        frame_Xsize = frame_Xsize - 1;
    end
    
    % To produce image data in final output form:
    % Loop through frames:
    for p = 1:numFrames
        % frame = imread(image_path,p); % read frame number p, read full image.
        frame = imread(image_path,p,'PixelRegion', {[1,frame_Ysize], [1,frame_Xsize]}); % read only a part of image, to have a square image at the end.
        frame = single(frame);  % to class single.
        % frame = double(frame);  % to class double.
        image_data(p).frame_data = frame; 
    end
       

end


%% For .mat files 

if isempty(matImageExists)==0
     
    load(image_path); % load mat file.
    % The mat file should contain a matrix of size frameSizeY x frameSizeX
    % x numFrames.
    
    file_info = whos('-file',image_path);
    % this returns a structure with fields, eg.:
    %     name: 'A_image_sequence'
    %     size: [400 400 100]
    %     bytes: 128000000
    %     class: 'double'
    %     global: 0
    %     sparse: 0
    %     complex: 0
    %     nesting: [1x1 struct]
    %     persistent: 0
    
    % Get useful info for final output:
    numFrames = file_info.size(3); % number of frames in sequence.
    frame_Ysize = file_info.size(1);
    frame_Xsize = file_info.size(2);
%     % Make the image be a square image:
%     frame_Ysize = min(frame_Ysize,frame_Xsize);
%     frame_Xsize = min(frame_Ysize,frame_Xsize);

    % If the frame size is an odd number:
    if mod(frame_Ysize,2)~=0 % modulus after division
        frame_Ysize = frame_Ysize - 1;
        frame_Xsize = frame_Xsize - 1;
    end
        
    % To produce image data in final output form:
    image_sequence = eval(file_info.name); % evaluate the string containing the matrix name to get the matrix data into a variable.
    % Loop through frames:
    for p = 1:numFrames
        frame = image_sequence(:,:,p);
        frame = single(frame);  % to class single.
        % frame = double(frame);  % to class double.
        image_data(p).frame_data = frame; 
    end
    
end



%% Output info to command window:

disp(['The total number of frames in this image sequence is: ',num2str(numFrames)]) %output total number of frames to command window.
disp(' ') % empty line
disp(['The chosen frame number is: ',num2str(frame_number)]) 
disp(' ') % empty line
disp(['The frame has size (YxX): ',num2str(frame_Ysize),' x ',num2str(frame_Xsize)])
disp(' ') % empty line

frame = image_data(frame_number).frame_data;
frame = double(frame);
imshow(frame,[],'Border','tight');

end