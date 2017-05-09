function [frame]=extract1frame(frame_number)
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
% This function chooses a .sif image sequence from an input dialog box, 
% returns the image data corresponding to its frame number "frame_number", and plots it.
% Note: only opens .sif images!
% See "extract_image_sequence_data.m" for a more general function to
% open/read image sequences.
% 
% example of how to call this function: extract1frame(100) extracts and plots
% frame 100 of the chosen image sequence and plots it.

% Data is in 'C:\Isabel\ExperimData\HeikoData\'
disp('Data is in C:\Isabel\ExperimData\HeikoData\')
disp(' ') % empty line

% uigetfile opens a file dialog box to choose image data file:
[file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
data_folder_path = strcat(path_data,file_data);
disp(data_folder_path) % write .sif image path to command window.

% first get size of .sif image file: (see IO_Input folder, SifFunctions.txt, or page 95 my notebook 1):
[ReturnCode, numFrames, ImageSize, TotalAcquisitionSize]=GetAndorSifSize(data_folder_path,0);
disp(' ') % empty line
disp(['The total number of frames in this image sequence is: ',num2str(numFrames)]) %output total number of frames to command window.
% numFrames is the length of the sequence, the number of frames.
% ImageSize is the size of the image, e.g. 512*512.
% TotalAcquisitionSize is numFrames*ImageSize

% read .sif image data. sifData is an array of 1x1 structures, as many
% columns as frames on image sequence. Reads frames 1 to numFrames:
[sifData] = read_sif_data_direct(data_folder_path,numFrames,ImageSize,1,numFrames);
% sifData is a cell with as many elements as image frames.


frame = sifData{frame_number}.sliceData; % selects the chosen frame number of the cell using {},
% and extracts frame data which is saved in the field 'sliceData'.
% % Note that the original frame up to here has intensity values which are larger than 1 and is of class single:
% min(min(frame))
% max(max(frame))
% class(frame)

% frame = mat2gray(frame);% transform matrix into a grayscale image (with values between 0 and 1) to do image operations on it.

frame = double(frame);
frame = imrotate(frame,90);
imshow(frame,[],'Border','tight');
