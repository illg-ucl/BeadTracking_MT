function plotBeadTrajNumbers(image_label,minPointsTraj) 
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
% This function is not finished, IN PROGRESS!!


%% Read in the image-sequence data:

% Read image-sequence file:
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.


%% Get path for trajectory data (excel file):

% You need to be in the correct directory before running the function!!!!
% Find paths in current folder which contain 'image_label' string:
trajXlsPath0 = dir(strcat('*',image_label,'*.xls')); % Trajectory data path (excel file with the full trajectories as returned by function "linkTrajSegments.m").
% Error control:
if isempty(trajXlsPath0) % If there is no .xls trajectory data file for such image number, show error and exit function:
    error('Check you are in the correct directory and run again. No .xls file found for that image number. Make sure image number is in between quotes ''.'); 
end
trajXlsPath = trajXlsPath0.name;
% trajXlsPath0 is a structure and the file names is stored in the field 'name'.


%% Open all trajectory data (excel file): 
% (.xls file previously generated with functions 'FindTrajects' and 'linkTrajSegments'):

% error control:
if 2~=exist('xlsread') 
    error('Check file dependencies - you need to install xlsread'); 
end

% Open excel file and read the data:
[NUMERIC,TXT,RAW]=xlsread(trajXlsPath,'Track results'); % import the data in the sheet named 'Track results'.
% Import the column heads and assign ID
colheads = TXT;
% The column titles are: CentreX, CentreY, ClipFlag, rsqFitX, rsqFitY,
% FrameNumber, BeadNumber,	TrajNumber.

% Generate ID: ID is a structure with fiels with the same names as the
% column titles, and each has an ID value of 1, 2, 3, etc (see below).
for i=1:numel(colheads) 
    ID.(colheads{i}) = find(strcmp(TXT,colheads{i})); 
end
% eg. ID = 
%         CentreX: 181.1559
%         CentreY: 393.0774
%        ClipFlag: 0
%         rsqFitX: 0.9997
%         rsqFitY: 0.9991
%     FrameNumber: 1
%      BeadNumber: 3
%      TrajNumber: 3

% The trajectory number column:
traj = NUMERIC(:,ID.TrajNumber); % NUMERIC is the numeric data read from the excel file (without the row of column titles).
disp('File Loaded successfully!');

% Get individual tracks:

% List the points at which the trajectory number first appears:
[A,I,J] = unique(traj,'first'); % [A,I,J] = UNIQUE(traj,'first') returns the vector I to index the first occurrence of each unique value in traj.  
% A has the same values as in traj but with no repetitions. A will also be sorted.
% List the points at which the trajectory number last appears:
[A,Y,Z] = unique(traj,'last'); % UNIQUE(traj,'last'), returns the vector Y to index the last occurrence of each unique value in traj.

% Get the number of tracks (no. of different trajectories):
numtracks = numel(A);

% Create tracks structure:
tracks(1:numtracks) = struct('trajNumber',[],'xvalues0',[],'yvalues0',[],'xvalues',[],'yvalues',[],'intensity',[],'I0_IspotFit',[],'sigma_IspotFit',[],'msd_unavg',[],'frame',[],'timeabs',[],'timerel',[],'numel',[],'minNumPointsInTraj',[],'deltaTime',[],'msd',[],'errorMsd',[],'errorMsdRelPercent',[],'disp',[],'SNR',[],'bg_noise_offset_afterBGsubtract',[],'BgNoiseStd',[],'IbgAvg',[],'IinnerTot',[],'rsqFit',[]);
del = []; % initialise for later.

for i=1:numtracks 
    % i
    a = I(i); % index for starting point in trajectory.
    b = Y(i); % index for ending point in trajectory.
    
    % Delete tracks that are less than minPointsTraj data points long:
    if b-a+1 >= minPointsTraj  % Only analyse tracks which have at least "minPointsTraj" points (frames) in them (5, or 15, e.g.).
    
    data{i} = NUMERIC(a:b,:);
    % tracks(i).XLS.track_index = A(i);
    tracks(i).trajNumber = A(i);
    % all values in pixels.
    tracks(i).xvalues0 = data{i}(1:end,ID.CentreX); % original xvalues in image (used later for plotting traj on image).
    tracks(i).yvalues0 = data{i}(1:end,ID.CentreY); % original xvalues in image (used later for plotting traj on image).
    % Set origin to zero:
    tracks(i).xvalues = tracks(i).xvalues0 - (tracks(i).xvalues0(1)); % xvalues relative to the first one in the trajectory.
    tracks(i).yvalues = tracks(i).yvalues0 - (tracks(i).yvalues0(1)); % % yvalues relative to the first one in the trajectory.
    tracks(i).msd_unavg = tracks(i).xvalues.^2+tracks(i).yvalues.^2; % squared displacement from the origin: x^2 + y^2.
    tracks(i).frame = data{i}(1:end,ID.FrameNumber); % frame number.
    tracks(i).timeabs = data{i}(1:end,ID.FrameNumber).*tsamp; % tsamp is the time between frames.
    tracks(i).timerel = tracks(i).timeabs-tracks(i).timeabs(1); % Set the first frame analysed as time zero reference (not used for now). 
    tracks(i).numel = b-a+1; % Number of points in the track. Isabel: it used to be b-a, I changed it to b-a+1.
    tracks(i).minNumPointsInTraj = minPointsTraj;
    tracks(i).rsqFitX = data{i}(1:end,ID.rsqFitX); % r-square of parabolic fit to centre peak of cross-correlation along x, for bead centre finding.
    tracks(i).rsqFitY = data{i}(1:end,ID.rsqFitY); % same along y.
    tracks(i) = getDisplacement(tracks(i),tsamp); % calculate msd and its error and add it to result structure.
    
    else
        % save indices to delete later:
        del(i) = i;     
    end
    
end

% Delete tracks which were too short: 
tracks(find(del))=[];

% ========================
%% For trajectory n n (finish!!)

n = 1; % for first frame:
% Show video of trajectory overlaid on actual image:
frames_list = analysedAllTraj(n).frame; % list of frame numbers in trajectory n.
x_values = analysedAllTraj(n).xvalues0; % list of original x centres of spots in trajectory n.
y_values = analysedAllTraj(n).yvalues0; % list of original y centres of spots in trajectory n.

%% For frame k:
k= 1; % use first frame:
frame = image_data(frames_list(k)).frame_data; % extract frame data which is stored in field 'frame_data'.
frame = double(frame);

imshow(frame,[],'Border','tight'); % show image scaled between its min and max values ([]).
hold on;

plot(x_values(k),y_values(k),'o','Color','g','MarkerSize',12) % plot accepted spot centres in green.
pause(0.1); % this pause is needed to give time for the plot to appear (0.1 to 0.3 default)
hold off;

% ==============

%% Plot frame average ROI with numbers of final accepted spots (index nn) overlaid 
% on top to be able to visually identify them:
figure; 
imshow(frame_Gray,[],'Border','tight');
hold on;
for nn=1:length(results)
    text(results(nn).Xcentre,results(nn).Ycentre,num2str(nn),'Color',[1 1 0],'FontSize',8); % number in yellow.
end
pause(0.1); % this pause is needed to give time for the plot to appear
hold off;
% SAVE current FIGURE as a .png and then close the figure window:
figName = strcat('Avg_frame_Accepted_Spot_Numbers',image_label); % choose name of figure file for saving .png.
saveFigurePNG(result_folder_name,figName); % See saveFigurePNG.m.
