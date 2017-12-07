function scriptToAnalyseManyVideos
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

% Move into correct data directory:
cd('C:\Users\isabel\ISABEL\UCL_Isabel\DataAnalysis\JohnnyNguyen_PSmagnetophoresis\Analysis PSBeads\x0')

%% PARAMETERS:

% Label used to produce output analysis excel files
data_set_label = 'analysis'; 

% Minimum number of points in track for it to be accepted:
minPointsTraj = 10;

% Make sure current directory is the directory containing the video files:
data_folder = cd;

% Choose video file extension to analyse:
videoFile_extension = '.avi';

% %% Regions to exclude from analysis:
% % Exclude certain region from all images in all videos.
% % Regions to exclude are given by start coordinates (x_start, y_start) for
% % top left corner of rectangle on image and end coordinates (x_end, y_end)
% % for bottom right corner of rectangle. Coordinates delimit rectangular
% % boxes on image and several rectangles can be defined. All rectangles are
% % put together in structure input excludedRegions:
% % Example, for Sonia's images we exclude the following regions:
% excludedRegions.list_xstart = [1 1 1 130];
% excludedRegions.list_xend = [820 112 100 215];
% excludedRegions.list_ystart = [581 492 1 1];
% excludedRegions.list_yend = [614 580 175 54];


%% Find video files with the chosen extension in current directory and make list of labels:

listVideoNames0 = dir(strcat('*',videoFile_extension));

% Error control:
if isempty(listVideoNames0) % If there are no video files of chosen extension, show error and exit function:
    error('Check you are in the correct directory. No video files of chosen extension found in folder.');
end

% List of video file names:
listVideoNames = {listVideoNames0.name}; % cell array of strings with video file names.

% Generate list of file labels, videoLabel, containing everything in the
% file name that comes before the file extension:
for k=1:length(listVideoNames)
    fullName = listVideoNames{k};
    pos = strfind(fullName,videoFile_extension); % position of the start of the string videoFile_extension in the file name.
    videoLabel{k} = fullName(1:(pos-1)); 
end



%% Loop through videos to analyse them:

for i = 1:length(videoLabel)   
   i
   videoLabel{i}
   % Track particles in videos:
   tracks{i} = FindTrajectsParticles(videoLabel{i},1,'end',excludedRegions);
   
   % save tracking results:
   save 'resultStructures.mat' 'tracks'
   
   linkTrajSegmentsParticles(videoLabel{i},1,'end',tracks{i},data_set_label);

   % Plot and save particle numbers on png:
   plotParticleTrajNumbers(videoLabel{i},data_set_label,10); 
   cd(data_folder);
   
   close all
end
    
 
%% Loop through videos to further analyse tracks:


for i = 1:length(videoLabel)   
   i
   videoLabel{i}
   % Go through tracks, decide which are "good" ones and save them to .mat
   % file in current directory for further processing:
   good_tracks{i} = goThroughParticleTracksVideo(videoLabel{i},data_set_label,1,'end',minPointsTraj,maxMajorAxisLength,0); 
   
   % Further process tracks to produce one plot and excel file per track:
   % showManyParticleTrajAnalysis(image_label,data_set_label,n_traj_start,n_traj_end,start_frame,tsamp,pixelsize_nm,showVideo,saveAvi,minPointsTraj)
   processedTrajs{i} = showManyParticleTrajAnalysis(videoLabel{i},data_set_label,1,'end',1,0.0333,1,1,1,minPointsTraj);
   % save processed results:
   save 'resultTracks.mat' 'processedTrajs'
end

