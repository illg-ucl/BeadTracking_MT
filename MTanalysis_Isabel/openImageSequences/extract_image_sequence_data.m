function [numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label)
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
%
% Read image sequence (.sif, .dv, .tif, .m4v or .mat data) and return image data in a useful
% form.
%
% INPUTS: 
% image_label: string that labels a given image sequence found in current
% folder. The code finds the path of the image file automatically based on a string label 
% that is equal to the file name without the file extension. For example,
% for image video file "210217r25.tif", image_label would be the
% string '210217r25'. Same throughout the entire BeadTracking_MT code.
%
% OUTPUTS:
% numFrames: number of frames in the image sequence.
% frame_Ysize: vertical (first dimension of matrix) size of image in pixels.
% frame_Xsize: horizontal (second dimension of matrix) size of image in pixels.
% image_data: structure array containing as many elements as frames. 
% Each element has field: frame_data: matrix data for each frame 
% with real intensity numbers (i.e., not between 0 and 1). 
% image_path: string containing path to image.
%
% To get frame number "p" do: image_data(p).frame_data.

%% Initial stuff

dvImagePath0 = dir(strcat(image_label,'.dv')); % Image sequence data path if the image is a .dv file.
sifImagePath0 = dir(strcat(image_label,'.sif')); % Image sequence data path if the image is a .sif file.
tifImagePath0 = dir(strcat(image_label,'.tif')); % Image sequence data path if the image is a .tif file.
m4vImagePath0 = dir(strcat(image_label,'.m4v')); % Image sequence data path if the image is a .m4v file.
matImagePath0 = dir(strcat(image_label,'.mat')); % Image sequence data path if the image is a .mat file.

% Error control:
% Sometimes we use a .mat file containing the good track numbers with a name which contains also
% the image sequence number and .mat, so in order to not think that is an
% image sequence, we include the following lines.
% If there is a .mat file and it contains the string "good_track_nums_" in
% its name, then:
if ~isempty(matImagePath0) && ~isempty(strfind(matImagePath0.name,'good_track_nums_'))
    matImagePath0 = [];
end

% Error control if neither a .sif, .tif, .mat or .dv file can be found with that label:
if isempty(dvImagePath0) && isempty(sifImagePath0) && isempty(tifImagePath0) && isempty(matImagePath0) % If there is no .sif, .tif or .dv image sequence file for such image_label, show error and exit function:
    error('Check you are in the correct directory and run again. No .sif, .tif, .dv or .mat file found for that image_label.');
end

% Turn off image size adjust warning: "Warning: Image is too big to fit on screen; displaying
% at X%". Warning identifier  is 'images:initSize:adjustingMag'.
warning('off','images:initSize:adjustingMag')


%% For .dv files 
% This needs bfopen.m and loci_tools.jar (see folder openDVfiles):

if isempty(dvImagePath0)==0
    
    image_path = dvImagePath0.name;

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
        % image data needs to be of class double for future operations (for findSpotCentre1frame.m).
    end
    
    % Add other useful info to final output:
    numFrames = series1_numFrames;
    frame_Ysize = size(image_data(1).frame_data,1);
    frame_Xsize = size(image_data(1).frame_data,2);

end


%% For .sif files:

if isempty(sifImagePath0)==0
    
    image_path = sifImagePath0.name;
    
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
        % image data needs to be of class double for future operations (for
        % findSpotCentre1frame.m).
    end
    
    % Get frame dimensions:
    frame_Ysize = size(image_data(1).frame_data,1);
    frame_Xsize = size(image_data(1).frame_data,2);
    
end


%% For .tif files 

if isempty(tifImagePath0)==0
    
    image_path = tifImagePath0.name;

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
        % image data needs to be of class double for future operations (for findSpotCentre1frame.m).
    end
    
    

end


%% For .m4v files 

if isempty(m4vImagePath0)==0
    
    image_path = m4vImagePath0.name;    
    m4v_info = VideoReader(image_path);
    % Add other useful info to final output:
    frame_Ysize = m4v_info.Height;
    frame_Xsize = m4v_info.Width;
%     % Make the image be a square image:
%     frame_Ysize = min(frame_Ysize,frame_Xsize);
%     frame_Xsize = min(frame_Ysize,frame_Xsize);

    % If by any chance the frame size is an odd number:
    if mod(frame_Ysize,2)~=0 % modulus after division
        frame_Ysize = frame_Ysize - 1;
        frame_Xsize = frame_Xsize - 1;
    end
    
    % To produce image data in final output form:
    % Read in frames:
    k = 1;
    while hasFrame(m4v_info)
        frame = readFrame(m4v_info);
        frame = single(frame);  % to class single.
        % frame = double(frame);  % to class double.
        % Convert RGB values to grayscale values by forming a weighted sum
        % of the R, G, and B components, see rgb2gray:
        image_data(k).frame_data = 0.2989 * frame(:,:,1) + 0.5870 * frame(:,:,2) + 0.1140 * frame(:,:,3);
        k = k+1;
    end
    
    numFrames = k-1; % number of frames in sequence.    
    
end


%% For .mat files 

if isempty(matImagePath0)==0
    
    image_path = matImagePath0.name;
    
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
        % image data needs to be of class double for future operations (for findSpotCentre1frame.m).
    end
    
end




%% Output info to command window:

% disp(' ') % empty line
disp(image_path) % write .sif image name to command window.
disp(['The total number of frames in this image sequence is: ',num2str(numFrames)]) %output total number of frames to command window.
disp(' ') % empty line

%  % To show a video:
%  % Loop through frames:
% for p = 1:numFrames
%     imshow(image_data(p).frame_data,[]);
%     pause(0.1);
% end


% -------------------
% Alternatively:
% uigetfile opens a file dialog box to choose data file:
% [file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
% strcat('data (.sif image):','  ',path_data,file_data)
% open a file dialog box to choose analysis file:
% [file_analysis,path_analysis] = uigetfile({'*.xls'}, 'Chose analysis file (trajectory):');
% strcat('analysis file (.xls trajectory):','
% ',path_analysis,file_analysis)