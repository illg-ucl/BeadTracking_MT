function processedTraj = showBeadTrajAnalysis(trajXlsPath,image_data,analysedAllTraj,n_traj,start_frame,tsamp,pixelsize_nm) 
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
% Analyse and process trajectory data (only traj number n_traj) for a given image sequence:
% - get plot and fit of exponential decay of Intensity vs time, 
% - get plot (and fit) of mean square displacement (msd) vs Delta t with error bars, 
% (- Interpolate Intensity vs time to have it uniformly sampled)
% - apply Chung-Kennedy filter to Intensity vs time data,
% - calculate pairwise differences of filtered (or unfiltered) data and
%   their distribution,
% - calculate power spectrum (Fourier transform) of pairwise difference
%   distribution,
% - find peaks in the power spectrum to get the intensity step.
%  
% NOTE: this function works for short tracks too. 
% It only fails if no. of points in track is < 3.
%
% The program asks the user to enter the parameter values for the
% Chung-Kennedy filter.
% Then show results only for the selected single trajectory (given by
% 'n_traj') and save that plot as a .png.
% 
% Inputs: 
% - 'trajXlsPath' is the path to an excel file with the full trajectories as returned by
% function "linkTrajSegmentsBeads.m".
% - 'image_data': structure array with each element being a 2D matrix with 
% the image data for each frame in the image sequence. 
% - 'analysedAllTraj': analysed trajectory information. It is the result
% from doing: analysedAllTraj = analyseTraj(trajXlsPath,tsamp,minPointsTraj); 
% It is a structure array with as many elements as
% analysed trajectories (the ones with at least 5 points in them), and
% with fields 'XLS', 'xvalues', 'yvalues', 'intensity', 'msd_unavg',
% 'timeabs', 'timerel', 'numel', 'deltaTime', 'msd', 'errorMsd',
% 'errorMsdRelPercent' and 'disp'.
% - 'n_traj' refers to the number of trajectory out of all the analysed
% trajectories for which the results are shown.
% - 'start_frame' is the number of frame considered as the origin of time (as t=0), in frames. It is the first frame for which
% the shutter is fully open and we detect fluorescence (written in my notebook when I analysed each image sequence).
% - 'tsamp' is the sampling time, used to calibrate the absolute time, to go
% from frames to time in seconds. 
% It is the time between frames in seconds. Use tsamp = 1 for the time to be in units of frames. A proper calibration
% have tsamp = 40*10^(-3), i.e., 40ms per frame, for example.
% start_frame*tsamp is therefore the absolute time origin in seconds.
% - 'pixelsize_nm': pixel size in nm (35.333nm for OXPHOS data).
%
% Output: 'processedTraj' is a cell array 4 elements:
% The first element, {1} contains summary results and is a structure with fields: 
%     'good_tracking_flag'
%     'XLSpath'
%     'minNumPointsInTraj'
%     'ImageSizeHorizPix'
%     'ImageSizeVerticPix'
%     'AnalysedTrajNum'
%     'OriginalTrajNum'
%     'FirstTrajFrame'
%     'LastTrajFrame'
%     'NumDataPoints'
%     'TrajStartTime'
%     'TrajEndTime'
%     'TrajDuration'
%     'TimeBetweenFrames'
%     'FrameForTimeOrigin'
%     'AbsTimeOrigin'
%     'pixelsize_nm'
%     'Track_meanX_0'
%     'Track_meanY_0'
%     'Track_meanX_cellOrigin'
%     'Track_meanY_cellOrigin'
%     'Track_meanXvalue_cellCoords'
%     'Track_meanYvalue_cellCoords'
%     'TopOrBottom'
%     'WindowWidthCKfilter'
%     'WeightExponentCKfilter'
%     'useFiltered'
%     'numBins'
%
% The next element, {2}, is a structure with the following fields, each containing a vector: 
%     'frame'
%     'timeabs'
%     'intensity'
%     'xvalues0'
%     'yvalues0'
%     'xvalues'
%     'yvalues'
%     'x_values_cell'
%     'y_values_cell'
%
% The next element, {3}, is a structure with the following fields, each containing a number: 
%     'I0_fit'
%     'stDev_I0'
%     'tau_fit'
%     'stDev_tau'
%     'rsq_fit_I'
%     'I0_fit_wo'
%     'stDev_I0_wo'
%     'Ioffset_fit_wo'
%     'stDev_Ioffset_wo'
%     'tau_fit_wo'
%     'stDev_tau_wo'
%     'rsq_fit_I_wo'
%     'I_line_offset'
%     'stDev_I_line_offset'
%     'I_line__slope'
%     'stDev_I_line__slope'
%     'I_line_rsq'
%     'Npoints_for_I_line_fit'
%
% The next element, {4}, is a structure with the following fields, each
% containing a vector: 
%     'deltaTime'
%     'msd'
%     'errorMsd'
%     'errorMsdRelPercent'
%
% The next element, {5}, is a structure with the following fields, each
% containing a vector: 
% 'BinCentresPos': x axis of histogram of Intensity pairwise differences.
% 'FreqCountsPos': y axis of histogram of Intensity pairwise differences.
%
% The next element, {6}, is a structure with the following fields, each
% containing a vector:
% 'IstepAxis' : x axis of spectrum (Fourier transform) of the previous histogram.
% 'PowerSpectrum' : y axis of spectrum (Fourier transform) of the previous histogram.
%
% The next element, {7}, is a structure with the following fields, each
% containing a vector: 
% 'Isteps':  vector with intensity steps, i.e., found peaks in spectrum.
% 'power_Fourier': vector with Fourier power (probability) of found
% peaks in spectrum.
% (see "showTrajAnalysis.m")
%
% The next element, {8}, is a structure with the following fields, each
% containing a number: 
%     'xsize_default_cell_region'
%     'ysize_default_cell_region'
%     'local_cell_region_xleft'
%     'local_cell_region_xright'
%     'local_cell_region_ytop'
%     'local_cell_region_ybottom'
%     'centre_of_mass_X'
%     'centre_of_mass_Y'
%     'coord_transform_angle'
%     'cell_width'
%     'cell_length'
%     'traj_in_pole'
%
% Note: function "showManyTrajAnalysis.m" calls this function "showTrajAnalysis.m".
%
% Further field explanation:
% 'XLSpath' % path of the excel file with the traj analysis (trajXlsPath)
% 'AnalysedTrajNum' % number of analysed trajectory.
% 'FirstTrajFrame' % First frame in this trajectory. 
% 'NumDataPoints' % number of data points in trajectory.
% 'FrameForTimeOrigin' % Frame corresponding to the time origin of the
%   image sequence (when shutter opens and camera starts recording
%   fluorescence).
% 'TimeBetweenFrames' % Time between frames in seconds
% 'AbsTimeOrigin' % Absolute time origin in seconds: it is the FrameForTimeOrigin times the
%   TimeBetweenFrames.
% 'IntensityAtTimeOrigin' % Integrated fluorescent spot intensity at the
%   time origin (from fit). 
% 'StDevI0' % error of the IntensityAtTimeOrigin (standard deviation of the
%   fit).
% 'Tau' % time constant of the exponential fit of the integrated spot
%   intensity versus absolute time.
% 'StDevTau' % error of Tau from the fit, as its standard deviation.
% 'rsqFitI' % r squared of exponential fit of the integrated spot
%   intensity versus absolute time.
% 'WindowWidthCKfilter' % width of window (in data points) for
%   Chung-Kennedy filter of integrated Intensity data.
% 'WeightExponentCKfilter' % Exponential factor for the weights in the same Chung-Kennedy filter. 
% 
%
% Use analysedAllTraj = analyseTraj(trajXlsPath,tsamp,minPointsTraj);
% and showTrajAnalysis(trajXlsPath,analysedAllTraj,n_traj,0,1) for quick results (time in frames, uncalibrated).
%
% E.g. of how to call this function: 
% pathTraj = 'C:\Isabel\ExperimData\HeikoData\ATPase-GFP\ATPase-GFP TIRF\ATPase-GFP_89fullTrajs.xls';
% sampling_time = 0.04; 
% t1 = analyseTraj(pathTraj,sampling_time,minPointsTraj);
% procTraj = showTrajAnalysis(pathTraj,t1,1,4,sampling_time)


%% Display some info first:

% Some needed track parameters-info:
firstFrameInTraj = analysedAllTraj(n_traj).frame(1); % first frame in trajectory.
lastFrameInTraj = analysedAllTraj(n_traj).frame(end); % last frame in trajectory.
traj_duration = tsamp*(lastFrameInTraj-firstFrameInTraj); % duration of track in seconds.
nPointsInTrack = analysedAllTraj(n_traj).numel;

list_spot_width = analysedAllTraj(n_traj).sigma_IspotFit; % vector containing fitted gaussian width to spot at all times in track.
meanSpotWidth = mean(list_spot_width); % Average spot width in pixels.

% t_origin is the absolute time origin in seconds:
t_origin = start_frame*tsamp;

disp(' ') % empty line
disp(['The analysed file is: ',trajXlsPath]) 
disp(['The number of long enough trajectories in the file is: ',num2str(length(analysedAllTraj))]) 
disp(['The time origin for this image sequence is frame no.: ',num2str(start_frame)]) 
disp(' ') % empty line
disp(['TRAJECTORY number: ',num2str(n_traj)]) 
disp(['The number of data points in this trajectory is: ',num2str(nPointsInTrack)]) 

%% Create new directory for saving trajectory-analysis results for this image sequence:

% Make new folder (new directory) to save trajectory analysis results:
pos1 = strfind(trajXlsPath,'fullTrajs.xls'); % position of the start of the string 'fullTraj.xls' in the xls input file name.
new_folder_name = trajXlsPath(1:(pos1-1)); % Take the name of the input excel file (with the end bit 'fullTraj.xls' removed) as the new folder name.
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(new_folder_name); % make new directory.


%% Read in the image-sequence data:

% Read image-sequence file: 
%[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.
% --------------------------------------------------------------


%% Calculate frame average:

numFrames = length(image_data); % number of frames in image sequence.
frame_Ysize = size(image_data(1).frame_data,1);
frame_Xsize = size(image_data(1).frame_data,2);

% Initialise frame accumulation in order to later calculate a frame average:
frame_accumul = zeros(frame_Ysize,frame_Xsize);

for k = start_frame:numFrames  % loop through frames.
    % Get frame data: 
    frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    % Accummulate frames to then calculate frame average:
    frame_accumul = frame_accumul + frame;
end

% Calculate frame average as the accumulation of all frames divided by the number of frames:
frame_avg = frame_accumul/(numFrames-start_frame+1);


%% PARAMETERS

% I have saved all parameter values in the file "paramsForShowTrajAnalysis2.m"
% in the current directory. Just calling the name of that file loads the
% parameter values into the workspace:
paramsForShowTrajAnalysis2
% In this way, we can save different parameter sets for different data sets
% in an easy manner and run two Matlabs at the same time working with different parameter sets.

% % PARAMETER r:
% % Number of frames to average over when calculating a frame average in
% % order to get a Signal Mask to distinguish cell region from background
% % region. Note it has to be at most equal to the total no. of frames minus
% % start_frame:
% % CHECK: check if only a few first frames (20) needed or more:
% r = numFrames-start_frame; % number of frames to average over, starting from start_frame. (20, numFrames-start_frame)
% 
% % CHECK:
% % PARAMETERs for default cell region around trajectory (for calculation of local cell coordinates):
% xsize_default_cell_region = 60; % horizontal size of region around cell. (70-140) 
% ysize_default_cell_region = 60; % horizontal size of region around cell. (70-140)
% 
% % PARAMETER: 
% % to find out initial intensity we fit only the first
% % "Npoints_for_I_line_fit" points of the Intensity vs time trace:
% Npoints_for_I_line_fit = 30; % (30).
% % Note that the min. number of points in the trajectory is given by
% % minPointsTraj (usually 15) in showManyTrajAnalysis.m.
% % If the no. of points in this traj is lower than the one above:
% if analysedAllTraj(n_traj).numel < Npoints_for_I_line_fit
%     Npoints_for_I_line_fit = analysedAllTraj(n_traj).numel;
% end
% 
% % PARAMETERs for Chung-Kennedy filter: 'Wfilter' and 'Rfilter' are the window
% % size and weighting exponent. See 'chungKennedyFilter.m'.
% % Wfilter = input('Enter width of fiter window in no. of points: '); % request user input
% % Rfilter = input('Enter filter weighting exponent: '); % request user input
% % or:
% Wfilter = 3; % (default 3-5) use for quick analysis
% Rfilter = 1; % use for quick analysis
% % Error control: 
% % Window size for Chung-Kennedy filter cannot be larger than number of data points minus one:
% Wfilter = min(Wfilter,(analysedAllTraj(n_traj).numel-1)); % choose the lowest.
%     
% % PARAMETER mobility_rsq_limit: value of rsquare of fit above which we
% % accept that the msd versus delta-time fit is either Brownian diffusion
% % (linear) or confined diffusion:
% mobility_rsq_limit = 0.4;
% % Maximum timeconstant in seconds from confined-trajectory fit (if time
% % constant is too large, the fit is actually linear...) for trajectory to
% % be labelled as confined diffusion:
% max_conf_timeconst = 100; 
% % Note: the guesses for the fits of the msd to a line or to a saturating
% % curve are given later and their values can make the fits fair or not.
% 
% % PARAMETER nbins: number of bins in histogram of intensity pair-wise differences.
% % See fullPwD.m.
% nbins = 500; % (50, 200)
% 
% % PARAMETER x_limit_spectr: max Intensity step size to plot as max value in
% % x-axis of power spectra for each spot...
% x_limit_spectr = 5000;
% 
% % PARAMETERS for flags:
% % Max number of frames in track for it to be flagged as "short" track:
% max_NumFramesForShortTrack = 10; 
% % Number of frames in track for it to be flagged as "long" track:
% min_NumFramesForLongTrack = 50; 
% % Number of frames in track for it to be flagged as "very long" track:
% min_NumFramesForVeryLongTrack = 120; 
% % Used for flaging trajectories which are good for extrapolating and
% % obtainin the initial intensity Istart at the first frame in the track.
% % Maximum no. of frames away from TimeOrigin frame of sequence:
% max_framesAwayFromTimeOrigin = 10;
% % For good exponential fit flags:
% min_rsq_forGoodExpFit = 0.7;
% max_tau_relativeError_forGoodExpFit = 30; % as percentage.
% 
% % Intensity level close to background level (we consider only ~12 steps
% % of photobleaching above background). Flag tracks with low enough
% % intensity levels above background to try and determine stoichiometry
% % later.
% % We flag up a track if it has points with intensity below these:
% lowEnoughI_limit_bottom = 10000; % for bottom channel, GFP, green.
% lowEnoughI_limit_top = 3000; % for top channel, mCherry, red.
 

%% Calculate Signal Mask to distinguish cell region from background region:
% Use frame average (of first frames only) to calculate signal mask. 

% Calculate average of first "r" frames:
% Initialise frame accumulation in order to later calculate a frame average:
frame_accumul_0 = zeros(frame_Ysize,frame_Xsize);

% r is the number of frames to average over, starting from start_frame.
% See PARAMETERS section.

for k = start_frame:start_frame+r  % loop through frames.
    % Get frame data: 
    frame = image_data(k).frame_data; % extract frame data, stored in the field 'frame_data'.
    frame = double(frame);
    % Accummulate frames to then calculate frame average:
    frame_accumul_0 = frame_accumul_0 + frame;
end

% Calculate frame average as the accumulation of all frames divided by the number of frames:
frame_avg_0 = frame_accumul_0/(r+1);

frame_avg_0_Gray = mat2gray(frame_avg_0); % The input to function "getCellMaskAndBoundary" needs to be a grayscale image:

% CHECK: use the full image (top and bottom) to threshold it and find
% signal mask. If you want to use a local threshold value within a selected
% cell region, comment the line out and uncomment the line were
% getCellMaskAndBoundary2 appears.
% Get SignalMask to know where cells are, to distinguish cells from background:
% [SignalMask CellBoundaryMask] = getCellMaskAndBoundary(frame_avg_0_Gray); 
% SignalMask is a matrix with 1 at positions where cells are and 0 at
% background.



%% Plot frame average with all trajectory points and select cell region for local cell coordinates:

% Display frame average and all trajectory positions overlaid on it:
% figure(200) % create figure number 200.
h = figure('Tag','cell_cords','position',[1000 100 800 800]); 
subplot(2,2,1);
imshow(frame_avg,[],'Border','tight'); % display frame average scaled between its min and max values ([]).
hold on;
% Trajectory data:
x_values = analysedAllTraj(n_traj).xvalues0; % list of original x centres of spots in trajectory n_traj.
y_values = analysedAllTraj(n_traj).yvalues0; % list of original y centres of spots in trajectory n_traj.

if mean(y_values) < (frame_Ysize/2)
    circleColour_string = 'r'; % plot circles in red if in top half of image.
else
    circleColour_string = 'g'; % plot circles in red if in bottom half of image.
end

for k = 1:length(x_values)
    plot(x_values(k),y_values(k),'.-','Color',circleColour_string,'MarkerSize',3) % plot track positions in green.
end
hold off;


% Default cell region:
xleft = round(mean(x_values) - xsize_default_cell_region/2);
xright = round(mean(x_values) + xsize_default_cell_region/2);
ytop = round(mean(y_values) - ysize_default_cell_region/2);
ybottom = round(mean(y_values) + ysize_default_cell_region/2);

% Error control:
% if box for default_cell_region too close to image edge, limit value to image size:
xleft = max(1,xleft); % this makes xleft >=1.
xright = min(xright,frame_Xsize); % this makes xright<=frame_Xsize. 
ytop = max(1,ytop); % this makes ytop >=1.
ybottom = min(ybottom,frame_Ysize); % this makes ybottom<=frame_Ysize. 

% Default cell region corners and centre (for calculation of local cell coordinates) for plot:
def_cell_region_x = [mean(x_values),xleft,xleft,xright,xright]; % centre and corners of default cell region, x values.
def_cell_region_y = [mean(y_values),ytop,ybottom,ytop,ybottom]; % centre and corners of default cell region, y values.

% CHECK: give option of user input or not.
% use_default_cell_region = input('Use plotted cell region for local cell coordinates (1 for yes, 0 for no)? :');
use_default_cell_region = 1;
local_cell_region = [xleft xright ytop ybottom];
% CHECK: comment out some of the following or not:
% Particular code for a given image with two cells, one left, one right:
if use_default_cell_region == 0
    % Uncomment one of the following options:
    %local_cell_region = input('Enter region around cell of interest as [xleft xright ytop ybottom] (for local cell coordinates):'); % request user input
    % --------------------------------
%     if mean(x_values) < 234 % left cells
%         if mean(y_values) < 317 % top cell
%             local_cell_region = [205 232 272 318];
%         else   % bottom cell:
%             local_cell_region = [204 228 316 367];
%         end
%     else % right cells
%         if mean(x_values) < 284 % cell on the left
%             local_cell_region = [238 283 288 308];
%         else   % cell on the right:
%             local_cell_region = [284 327 270 301];
%         end
%     end
    % --------------------------------
    %     if mean(x_values) > 330 % right cell
    %         local_cell_region = [370 417 325 383]; % values entered by hand for a given image sequence.
    %     else % cells on the left
    %         if mean(y_values) < 326 % cell on the top
    %             local_cell_region = [208 245 277 326];
    %         else   % cell on the right:
    %             if mean(y_values) < 387 % cell in the middle
    %                 local_cell_region = [224 245 332 379];
    %             else % cell in the bottom
    %                 local_cell_region = [196 247 397 451];
    %             end
    %         end
    %     end
    % --------------------------------
    % local_cell_region = [205 264 152 209];
    % --------------------------------
%     if mean(x_values) < 250 % cell on the left
%         local_cell_region = [192 244 149 200];
%     else % cell on the right
%         local_cell_region = [289 349 17 56];
%     end
    % --------------------------------
    if mean(x_values) < 68
        local_cell_region = [21 73 62 105];
    else % top cells
        if mean(x_values) < 150
            local_cell_region = [70 113 116 149];
        else
            local_cell_region = [200 251 139 202];
        end
    end
    % --------------------------------
%     if mean(y_values) < 256 % red, top channel
%         if mean(x_values) < 180
%             if mean(y_values) < 108
%                 local_cell_region = [76 116 78 107];
%             else
%                 local_cell_region = [75 126 107 143];
%             end
%         else
%             if mean(y_values) < 78
%                 local_cell_region = [273 304 43 78];
%             else
%                 local_cell_region = [251 300 79 120];
%             end
%         end
%     else % bottom, green channel
%         if mean(x_values) < 180
%             if mean(y_values) < 370
%                 local_cell_region = [68 106 338 370];
%             else
%                 local_cell_region = [59 116 370 408];
%             end
%         else
%             if mean(y_values) < 338
%                 local_cell_region = [267 304 298 338];
%             else
%                 local_cell_region = [253 391 338 377];
%             end
%         end
%     end
    % --------------------------------
    xleft = local_cell_region(1);
    xright = local_cell_region(2);
    ytop = local_cell_region(3);
    ybottom = local_cell_region(4);
end


% CHECK:
% Use a local threshold within the selected cell region to threshold the
% image and find the cell-signal mask: use getCellMaskAndBoundary2.m instead
% of getCellMaskAndBoundary.m.
[SignalMask CellBoundaryMask] = getCellMaskAndBoundary2(frame_avg_0_Gray,local_cell_region); 
% frame_avg_0_Gray is the complete greyscale image.


subplot(2,2,3);
imshow(SignalMask,[],'Border','tight'); % display signal mask around cell region.
hold on;
% Plot default cell region corners and centre (for calculation of local
% cell coordinates):
plot(def_cell_region_x,def_cell_region_y,'+r','MarkerSize',3) 
% Plot new cell region corners in green (for calculation of local cell coordinates):
def_cell_region_x = [xleft,xleft,xright,xright]; % corners of default cell region, x values.
def_cell_region_y = [ytop,ybottom,ytop,ybottom]; % corners of default cell region, y values.
plot(def_cell_region_x,def_cell_region_y,'+g','MarkerSize',3) 
hold off;



%% Get local cell coordinates: 

% Coordinate transformation from (x,y) to (x', y'):
% Use cell mask and selected local cell region to find out coordinate transformation to go to local cell coordinates:
% Centre of mass:
% Xpos is a matrix of the same size as frame, containing x values for all
% pixels and similarly for Ypos. For the whole image:
[Xpos,Ypos] = meshgrid(1:frame_Xsize,1:frame_Ysize);
% SignalMask_skeleton = bwmorph(SignalMask,'skel',Inf); % apply operation skeleton infinite times til we just get line skeleton.
% SignalMask_skeleton = bwmorph(SignalMask,'skel',5); % apply operation skeleton 5 times. 
% CHECK: check size of disk is good for type of data (depends on size of cells) and check that method works well:
se = strel('disk',7);
SignalMask_eroded = imerode(SignalMask,se); % erode image with a disk of size 9.
SignalMask_skeleton = bwmorph(SignalMask_eroded,'skel',Inf); % get skeleton (line) of eroded image.
% SignalMask_skeleton = bwmorph(SignalMask,'skel',Inf); % get skeleton (line) of signal mask image.

Cx = SignalMask.*Xpos; % for whole signal mask.
Cy = SignalMask.*Ypos;
Cx_1 = Cx(ytop:ybottom,xleft:xright); % use cell region only.
Cy_1 = Cy(ytop:ybottom,xleft:xright); % use cell region only.
Cx_2 = Cx_1(Cx_1~=0); % extract non-zero values. X values of cell region for selected region of interest.
Cy_2 = Cy_1(Cy_1~=0); % extract non-zero values. Y values of cell region for selected region of interest.
% Centre of mass of selected cell region:
centre_of_mass_X = mean(Cx_2);
centre_of_mass_Y = mean(Cy_2);

Cx_skel = SignalMask_skeleton.*Xpos; % for skeleton line of signal mask only
Cy_skel = SignalMask_skeleton.*Ypos;
Cx_1_skel = Cx_skel(ytop:ybottom,xleft:xright); % use cell region only.
Cy_1_skel = Cy_skel(ytop:ybottom,xleft:xright); % use cell region only.
Cx_2_skel = Cx_1_skel(Cx_1_skel~=0); % extract non-zero values. X values of cell region for selected region of interest.
Cy_2_skel = Cy_1_skel(Cy_1_skel~=0); % extract non-zero values. Y values of cell region for selected region of interest.

% Change origin of coordinate system to the centre of mass of the cell
% region:
x_values_2 = x_values-centre_of_mass_X; % trajectory data is x_values, y_values.
y_values_2 = y_values-centre_of_mass_Y; % column vectors.
% Transform also coordinates of the cell region mask:
Cx_3 = Cx_2 - centre_of_mass_X;
Cy_3 = Cy_2 - centre_of_mass_Y;
Cx_3_skel = Cx_2_skel - centre_of_mass_X;
Cy_3_skel = Cy_2_skel - centre_of_mass_Y;

% angle = atan(Cy_3./Cx_3); % angle (deg) for each point in the cell mask.
% radial = sqrt(Cx_3.^2+Cy_3.^2);

% Fit Y coordinates versus X coordinates (of previous cell mask region, of its skeletons) to a line to find angle for coordinate transformation:
linearFunction = fittype('A + B*x','independent','x'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');
warning('off','curvefit:fit:noStartPoint');
% Guesses for fit parameters:
% guess_A = 0; 
% guess_B = 0; 
% options.StartPoint = [guess_A guess_B]; % give guess parameters for fit. 
% Use coeffnames(fit_result_I) later to find out order of parameters.
% Use warning('query','last') on command window to find warning identifier.

try % error control, catch error
    [fit_result gof] = fit(Cx_3_skel,Cy_3_skel,linearFunction,options);
    fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'A', second one is 'B'.
    % Use coeffnames(fit_result) to find out order of parameters.
    cell_region_offset = fit_param_values(1); % offset from fit.
    cell_region_slope = fit_param_values(2); % slope from fit.
catch ME1 % if there was an error in previous command:
    cell_region_offset = 0; 
    cell_region_slope = 0; 
end
% the slope is the tan of the fit angle. The angle for the transformation is minus the angle from the fit:
coord_transform_angle = -atan(cell_region_slope); % the angle is minus the arctan of the slope, in radians here.
% Coordinate transform, rotation matrix, rotation by angle coord_transform_angle:
rotation_matrix = [cos(coord_transform_angle) -sin(coord_transform_angle); sin(coord_transform_angle) cos(coord_transform_angle)]; % coordinate transformation matrix, rotation by angle coord_transform_angle.
% Calculate trajectory in new rotated coordinate system, local cell coordinates:
x_values_cell = zeros(length(x_values_2),1); % initialise size, column vector. 
y_values_cell = zeros(length(y_values_2),1); % initialise size, column vector.
for i=1:length(x_values_2) % make coordinate transformation, rotation:
    new_coords = rotation_matrix*[x_values_2(i); y_values_2(i)]; % new coordinates.
    x_values_cell(i) = new_coords(1); % trajectory-x in new cell-coordinate system (cell's long axis).
    y_values_cell(i) = new_coords(2); % trajectory-y in new cell-coordinate system (cell's transverse axis).
end
% transform also the coordinates of the cell region mask:
for i=1:length(Cx_3) % make coordinate transformation, rotation:
    new_coords = rotation_matrix*[Cx_3(i); Cy_3(i)]; % new coordinates.
    Cx_4(i) = new_coords(1);
    Cy_4(i) = new_coords(2);
end

% plot original coordinates of selected cell region: Y coordinates versus X
% coordinates:
subplot(2,2,2);
% On image plane original x coord grows from left to right and y coord
% grows from top to bottom, so for cell to look the same as in image, we
% plot -y coord:
plot(Cx_3,-Cy_3,'*b','MarkerSize',2) % plot positions in cell region mask as blue stars.
hold on;
plot(0,0,'+r','MarkerSize',8) % plot centre point.
% Cx_3,-Cy_3 coords have their origin in the centre of mass.
% Plot linear fit:
plot(Cx_3,-(cell_region_offset + cell_region_slope.*Cx_3),'k-') % black line
set(gca,'DataAspectRatio',[1 1 1])
ylim([1.3*min(Cy_3) 1.3*max(Cy_3)])
hold off; 

subplot(2,2,4);
plot(x_values_cell,-y_values_cell,'.-g','MarkerSize',4) % plot new trajectory coordinates as green circles.
hold on;
plot(Cx_4,-Cy_4,'b.','MarkerSize',4) % plot new coordinates of cell region mask.
set(gca,'DataAspectRatio',[1 1 1])
hold off;


% In the new local cell coordinate system, the x'-axis is the long cell
% axis, and the y'-axis is the transverse short cell axis.
% The cell centre is at the origin of the new cell-coordinate system.
% cell width:
cell_width = 2*mean([abs(min(Cy_4)) abs(max(Cy_4))]);
% cell length:
cell_length = 2*mean([abs(min(Cx_4)) abs(max(Cx_4))]);
% Pole or not:
% We consider the trajectory is in a cell-pole region if it has an average horizontal cell-coordinate (x') fulfilling:
if abs(mean(x_values_cell))>(cell_length/2-cell_width/2) % if mean value of trajectory y' position is close to edge of cell
    traj_in_pole = 1;
else
    traj_in_pole = 0;
end

figure(h) % bring figure to front on screen.

% CHECK: if user input is required or not:
disp(' '); % empty line.
disp(['traj_in_pole (1 for pole, 0 for no pole) is now: ',num2str(traj_in_pole)]);
% agree_trajInPole = input('Do you agree with traj_in_pole? (1 for yes, 0 for no): '); % request user input
agree_trajInPole = 1;
if agree_trajInPole == 0 % if user does not agree with pole/no-pole:
    traj_in_pole = not(traj_in_pole); % take opposite;
end

% Save results related to local cell coordinate system (all in pixels):
results_cell_coords.xsize_default_cell_region = xsize_default_cell_region;
results_cell_coords.ysize_default_cell_region = ysize_default_cell_region;
results_cell_coords.local_cell_region_xleft = xleft;
results_cell_coords.local_cell_region_xright = xright;
results_cell_coords.local_cell_region_ytop = ytop;
results_cell_coords.local_cell_region_ybottom = ybottom;
results_cell_coords.centre_of_mass_X = centre_of_mass_X;
results_cell_coords.centre_of_mass_Y = centre_of_mass_Y;
results_cell_coords.coord_transform_angle = coord_transform_angle; % note that the angle in the image plane has oposite sign to the angle in an x-y plot, because image plane has reflected y coord.
results_cell_coords.cell_width = cell_width;
results_cell_coords.cell_length = cell_length;
results_cell_coords.cell_width_nm = pixelsize_nm*cell_width;
results_cell_coords.cell_length_nm = pixelsize_nm*cell_length;
results_cell_coords.traj_in_pole = traj_in_pole;
% Trajectory in cell-coordinates is a vector, so it is saved in a different
% excel sheet.

% close(findobj('Tag','cell_cords')); % close figure 'cell_cords'.

%% Save graphical analysis results for local cell coordinate system (last figure):

% Move into folder previously created to save traj analysis results and
% SAVE current FIGURE (results for local cell coordinate system) as a .png and then close the figure window:
figName = strcat(new_folder_name,'_cellCoords_traj',num2str(n_traj)); % choose name of figure file for saving .png.
saveFigurePNG(new_folder_name,figName); % See saveFigurePNG.m.


%% Fit intensity data to an exponentially decaying function (with no offset):

% We fit IspTot, integrated spot intensity after background subtraction.
IforFit = analysedAllTraj(n_traj).intensity; % Intensity is assumed to be already background corrected (see IspTot in 'findSpotCentre1frame.m').
tforFit = analysedAllTraj(n_traj).timeabs; % absolute time (comes from frame number and time between frames).
exp_no_offset = fittype('I0*exp(-(t-tstart)/tau)','independent','t','problem','tstart'); % define exponential funtion to fit to, with 't' as independent variable and 'tstart' as a fixed parameter (constant);
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_I0 = mean(IforFit(1:min(analysedAllTraj(n_traj).numel,5))); % Estimate as mean value of first five points (20*10^4 for oxphos data).
guess_tau = 40*tsamp; % around 40 frames
% Use coeffnames(fit_result_I) later to find out order of parameters.
options.StartPoint = [guess_I0 guess_tau]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [0 0]; % Lower bounds for fit parameters. In order: I0, Ioffset, tau.
options.Upper = [Inf 30]; % Upper bounds for fit parameters. In order: I0, Ioffset, tau.
[fit_result_I gof] = fit(tforFit,IforFit,exp_no_offset,options,'problem',t_origin); % fit_result_I contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
% 't_origin' is the absolute time origin in seconds.
% fit_param_names = coeffnames(fit_result_I); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
fit_param_values = coeffvalues(fit_result_I); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
I0_fit = fit_param_values(1); % I0 intensity value from fit.
tau_fit = fit_param_values(2); % tau from fit.
rsq_fit_I = gof.rsquare; % rsquare coefficient of fit.
errors = confint(fit_result_I,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
stDev_I0 = errorSTDEV(1);
stDev_tau = errorSTDEV(2);

disp(' ') % empty line
disp('Intensity vs time exponential fit (with no offset) result: ') 
disp([' I0 = ',num2str(I0_fit),' +- ',num2str(stDev_I0),';   tau = ',num2str(tau_fit),' +- ',num2str(stDev_tau),' s.']) 

% results from exponential fit with no offset of intensity vs time:
results_I_fits.I0_fit = I0_fit;
results_I_fits.stDev_I0 = stDev_I0;
results_I_fits.I0_fit_percentError = 100*stDev_I0/I0_fit;
results_I_fits.tau_fit = tau_fit;
results_I_fits.stDev_tau = stDev_tau;
results_I_fits.tau_fit_percentError = 100*stDev_tau/tau_fit;
results_I_fits.rsq_fit_I = rsq_fit_I;


%% Fit intensity data to an exponentially decaying function (with offset "_wo"):

% IforFit = analysedAllTraj(n_traj).intensity; % Intensity is assumed to be already background corrected (see IspTot in 'findSpotCentre1frame.m').
% tforFit = analysedAllTraj(n_traj).timeabs; % absolute time (comes from frame number and time between frames).
exp_with_offset = fittype('I0*exp(-(t-tstart)/tau)+Ioffset','independent','t','problem','tstart'); % define exponential funtion to fit to, with 't' as independent variable and 'tstart' as a fixed parameter (constant);
options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% Guesses for fit parameters:
guess_I0 = mean(IforFit(1:min(analysedAllTraj(n_traj).numel,5))); % Estimate as mean value of first five points (20*10^4 for oxphos data).
guess_Ioffset = 0;
guess_tau = 40*tsamp; % around 40 frames
% Use coeffnames(fit_result_I) later to find out order of parameters.
options.StartPoint = [guess_I0 guess_Ioffset guess_tau]; % give guess parameters for fit. This avoids a warning message. Give in right order!.
options.Lower = [0 -1000 0]; % Lower bounds for fit parameters. In order: I0, Ioffset, tau.
options.Upper = [Inf 30000 30]; % Upper bounds for fit parameters. In order: I0, Ioffset, tau.
try % error control in case fit fails
    [fit_result_I_wo gof] = fit(tforFit,IforFit,exp_with_offset,options,'problem',t_origin); % fit_result_I_wo contains the fit coefficient values and their confidence intervals and "gof" gives the "good of fitness".
    % fit_param_names = coeffnames(fit_result_I_wo); % fit parameter names: needed to check once their order: first one is 'I0', second one is 'tau'.
    fit_param_values = coeffvalues(fit_result_I_wo); % parameter values resulting from fit. First one is 'I0', second one is 'tau'.
    I0_fit_wo = fit_param_values(1); % I0 intensity value from fit.
    Ioffset_fit_wo = fit_param_values(2); % Ioffset from fit.
    tau_fit_wo = fit_param_values(3); % tau from fit.
    rsq_fit_I_wo = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_result_I_wo,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    stDev_I0_wo = errorSTDEV(1);
    stDev_Ioffset_wo = errorSTDEV(2);
    stDev_tau_wo = errorSTDEV(3);
catch ME1
    fit_result_I_wo = [0];
    I0_fit_wo = [];
    Ioffset_fit_wo = [];
    tau_fit_wo = [];
    rsq_fit_I_wo = [];
    errors = [];
    errorSTDEV = [];
    stDev_I0_wo = [];
    stDev_Ioffset_wo = [];
    stDev_tau_wo = [];
end


disp(' ') % empty line
disp('Intensity vs time exponential fit (with offset) result: ') 
disp([' I0_wo = ',num2str(I0_fit_wo),' +- ',num2str(stDev_I0_wo),';   tau_wo = ',num2str(tau_fit_wo),' +- ',num2str(stDev_tau_wo),' s.',';   Ioffset_wo = ',num2str(Ioffset_fit_wo),' +- ',num2str(stDev_Ioffset_wo)]) 

% results from exponential fit with offset of intensity vs time:
results_I_fits.I0_fit_wo = I0_fit_wo;
results_I_fits.stDev_I0_wo = stDev_I0_wo;
results_I_fits.I0_fit_percentError_wo = 100*stDev_I0_wo/I0_fit_wo;
results_I_fits.Ioffset_fit_wo = Ioffset_fit_wo;
results_I_fits.stDev_Ioffset_wo = stDev_Ioffset_wo;
results_I_fits.Ioffset_fit_percentError_wo = 100*stDev_Ioffset_wo/Ioffset_fit_wo;
results_I_fits.tau_fit_wo = tau_fit_wo;
results_I_fits.stDev_tau_wo = stDev_tau_wo;
results_I_fits.tau_fit_percentError_wo = 100*stDev_tau_wo/tau_fit_wo;
results_I_fits.rsq_fit_I_wo = rsq_fit_I_wo;
results_I_fits.Ioffset_relativeTo_I0_percent = 100*Ioffset_fit_wo/I0_fit_wo;


%% Fit first points in intensity data to a line to get initial intensity:

% Fit only the first Npoints_for_I_line_fit points (see PARAMETERS section).
% IforFit = analysedAllTraj(n_traj).intensity; % Intensity is assumed to be already background corrected (see IspTot in 'findSpotCentre1frame.m').
% tforFit = analysedAllTraj(n_traj).timeabs; % absolute time (comes from frame number and time between frames).
IforFit_1 = IforFit(1:Npoints_for_I_line_fit);
tforFit_1 = tforFit(1:Npoints_for_I_line_fit);
% The initial intensity is the one at t=t_origin (tstart), i.e., the offset of the fit, 'A':
linearFunction2 = fittype('A + B*(t-tstart)','independent','t','problem','tstart'); % define linear funtion to fit to, with 't' as independent variable and 'tstart' as a fixed parameter (constant);
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
guess_A = mean(IforFit_1(1:min(analysedAllTraj(n_traj).numel,5)));
guess_B = 0; 
options.StartPoint = [guess_A guess_B];
[fit_result_I_line gof] = fit(tforFit_1,IforFit_1,linearFunction2,options,'problem',t_origin); 
fit_param_values = coeffvalues(fit_result_I_line); % parameter values resulting from fit. First one is 'A', second one is 'B'. 
% Use coeffnames(fit_result_I_line) to find out order of parameters.
I_line_offset = fit_param_values(1); % offset from fit. This is the initial intensity, I at t_origin or tstart.
I_line_slope = fit_param_values(2); % slope from fit.
I_line_rsq = gof.rsquare; % rsquare coefficient of fit.
errors = confint(fit_result_I_line,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
stDev_I_line_offset = errorSTDEV(1);
stDev_I_line_slope = errorSTDEV(2);

disp(' ') % empty line
disp('Fit first points in intensity data to a line to get initial intensity: ') 
disp([' I_line_offset = ',num2str(I_line_offset),' +- ',num2str(stDev_I_line_offset),';   I_line_slope = ',num2str(I_line_slope),' +- ',num2str(stDev_I_line_slope)]) 

% results from linear fit to first intensity points:
results_I_fits.I_line_offset = I_line_offset; % This is the initial intensity, I at t_origin or tstart.
results_I_fits.stDev_I_line_offset = stDev_I_line_offset;
results_I_fits.I_line_slope = I_line_slope;
results_I_fits.stDev_I_line_slope = stDev_I_line_slope;
results_I_fits.I_line_rsq = I_line_rsq;
results_I_fits.Npoints_for_I_line_fit = Npoints_for_I_line_fit;


%% Interpolate intensity data to have it uniformly sampled:
% This is not really used.

% % Interpolate intensity trace to have it uniformly sampled in time for
% % further processing (since traj data can be non-uniformly sampled):
% x1 = analysedAllTraj(n_traj).frame; % column vector of original frame numbers (possibly non-uniformly sampled).
% y1 = IforFit; % column vector of original intensities (possibly non-uniformly sampled).
% % IforFit is the corresponding vector of IspTot values, the integrated,
% % background-corrected intensities.
% x1min = min(x1);
% x1max = max(x1);
% x2 = (x1min:1:x1max)'; % uniformly sampled frame numbers, one by one.
% y2 = interp1q(x1,y1,x2); % interpolated intensity: one point for each value of x2.
% 
% % Plot interpolated intensity trace as an aid for the user to choose the
% % appropriate Chung-Kennedy filter parameters (next step):
% figure(100); % figure number 100, far enough from the usual ones produced: 1, 2, 3...
% plot(x1,y1,'g-*'); 
% hold on;
% plot(x2,y2,'bo');
% xlabel('frame number');
% ylabel('Intensity')
% 
% 
% frames_afterInterp = x2; % frame numbers evenly sampled.
% times_afterInterp = x2*tsamp; % times evenly sampled (in seconds).
% IspTot_afterInterp = y2; % IspTot values after interpolation, uniformly sampled.
% 
% % interpolate_data = input('use interpolated intensity trace instead of original data? (enter 1 if yes) '); % request user input
% % disp(' ') % empty line

% ------------------------------
% Let us not interpolate:
interpolate_data = 0;
% ------------------------------

if interpolate_data == 1 
    I_to_use = IspTot_afterInterp; 
    t_to_use = times_afterInterp;
    frames_to_use = frames_afterInterp;
else
    I_to_use = IforFit;
    t_to_use = tforFit;
    frames_to_use = analysedAllTraj(n_traj).frame;
end


%% Apply Chung-Kennedy filter to intensity trace (IspTot_afterInterp):

% 'Wfilter' and 'Rfilter' are the window size and weighting exponent of the
% Chung-Kennedy filter. See 'chungKennedyFilter.m'.

% Go to PARAMETERS section to find values of Wfilter and Rfilter.

[filtered_I,tx,dx,sd,dsd,xpre] = chungKennedyFilter(I_to_use,Wfilter,Rfilter); 
% window size of Wfilter(3) points and weighting factor of Rfilter(1). 
% Output 'filtered_I' is the filtered vector (ignore other outputs).

% % continue with figure(1), aid plot shows filtered trace to decide if we
% % use filtered data or original data:
% plot(frames_to_use,filtered_I,'black-','LineWidth',2)
% hold off;


%% Choose between filtered (or unfiltered) intensity data to calculate pair-wise differences:

% use_filtered = input('use filtered trace instead of original data? (enter 1 if yes) '); % request user input
use_filtered = 0; % use for quick analysis
% disp(' ') % empty line
% close(100) % close window of aid figure of interpolated intensity trace.
% continues below in part 2.


%% Start multifigure here to plot/display all results:

% Create figure. 'position' vector is [left, bottom, width, height]. 
h1 = figure('Tag','Trajectory results','color','white','units','inches','position',[4 2.7 16 9]); 
% left, bottom control the position at which the window appears when it
% pops. 

% Plot xvalues and yvalues (origin is first point in trajectory) versus time:
subplot(3,3,1);
plot(analysedAllTraj(n_traj).timeabs,analysedAllTraj(n_traj).xvalues,'.-b'); 
hold on;
plot(analysedAllTraj(n_traj).timeabs,analysedAllTraj(n_traj).yvalues,'.-g');
xlabel('t_{abs} (s)'); 
ylabel('x-x0 (b),  y-y0 (g) (pix)');
xlim([0 max(analysedAllTraj(n_traj).timeabs)]);
% legend('x-x0 (pix)','y-y0 (pix)');
legend('hide');
% display image file path:
axes_limits = axis; % get axes limits: axes_limits = [xmin xmax ymin ymax].
text(axes_limits(1),1.2*axes_limits(4),trajXlsPath,'Interpreter','none'); % display path of image sequence at top left of figure at position (xmin,ymax) of axes.
hold off;

% Plot trajectory in x-y plane (origin is first point in trajectory), axes in nm now:
subplot(3,3,4);
plot(pixelsize_nm.*analysedAllTraj(n_traj).xvalues,pixelsize_nm.*analysedAllTraj(n_traj).yvalues,'.-k');
xlabel('x-x0 (nm)'); 
ylabel('y-y0 (nm)'); 

% Intensity vs time plot:
% Original IspTot intensity values versus time (before interpolating):
subplot(3,3,2);
plot(tforFit,IforFit,'gX'); 
hold on; 
plot(fit_result_I,'b'); % plot exponential fit (no offset) as blue line.
plot(fit_result_I_wo,'r'); % plot exponential fit (with offset) as red line.
plot(fit_result_I_line,'c'); % plot linear fit to first points as light blue (cyan) line.
% Plot of filtered intensity: after applying Chung-Kennedy filter on the interpolated intensity trace:
if ~isempty(filtered_I)
    plot(t_to_use,filtered_I,'black');
end
hold off; 
% legend('data','exp fit, no offset','exp fit, with offset','linear fit 1st-points','filtered data');
legend('hide');
xlabel('t_{abs} (s)'); ylabel('Intensity (arb)'); xlim([0 max(analysedAllTraj(n_traj).timeabs)]); ylim([0 Inf]);  


%% Msd (mean square displacement) versus delta-time:

% Original msd data:
xdata_msd_0 = double(analysedAllTraj(n_traj).deltaTime); % delta-time is in seconds.
ydata_msd_0 = double(analysedAllTraj(n_traj).msd); % fits fail later if data is not of class double. msd is in pixels^2 here.
% Error bars for msd versus Delta t plot:
halfErrorBars_0 = double(analysedAllTraj(n_traj).errorMsd); % These are the lower and upper bounds of the error bars around the msd result points for the graph, the standard deviation.
msd_relative_errors_0 = double(analysedAllTraj(n_traj).errorMsdRelPercent);
% % --------------
% % Plot all Msd data:
% subplot(3,3,5);
% errorbar(xdata_msd_0,ydata_msd_0,halfErrorBars,halfErrorBars,'.r'); 
% xlabel('\Deltat (s)'); 
% ylabel('msd (pix^2)'); 
% xlim([0 max(xdata_msd_0)]); 
% ylim([0 1.5*max(ydata_msd_0)]);  
% hold on;

% Only fit data points with relative error < 150%:
% % Use warning('query','last') to find out error identifier:
% warning('off','MATLAB:ignoreImagPart'); % turn off warning message for Warning: "Imaginary part ignored in comparison operation". 
% Select only real numbers (for large errors, msd_relative_errors is Inf and values are imaginary...)
pos_acceptedError_0 = []; % initialise empty vector.
for i=1:length(msd_relative_errors_0)
    if isreal(msd_relative_errors_0(i))==1        
        pos_acceptedError_0 = [pos_acceptedError_0  i]; % positions of real numbers in vector.
    end
end
msd_relative_errors_bis = msd_relative_errors_0(pos_acceptedError_0); % vector of relative errors excluding Inf and imaginary numbers. 
xdata_msd_bis = xdata_msd_0(pos_acceptedError_0);
ydata_msd_bis = ydata_msd_0(pos_acceptedError_0);
halfErrorBars_bis = halfErrorBars_0(pos_acceptedError_0);
% Now take only values with relative error < 150% :
pos_acceptedError = find(msd_relative_errors_bis < 150);
xdata_msd = xdata_msd_bis(pos_acceptedError);
ydata_msd = ydata_msd_bis(pos_acceptedError);
halfErrorBars = halfErrorBars_bis(pos_acceptedError);
msd_relative_errors = msd_relative_errors_bis(pos_acceptedError);

% -----------------------
% Fit the MSD versus delta-time to a LINE (Brownian-linear mobility):
linearFunction = fittype('A + B*x','independent','x'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% CHECK: values of guesses:
guess_A = 0; % offset guess in pix^2.
guess_B = 5; % slope guess in pix^2/s.
options.StartPoint = [guess_A guess_B];
try % error control
    [fit_msd_line gof] = fit(xdata_msd,ydata_msd,linearFunction,options);
    fit_param_values = coeffvalues(fit_msd_line); % parameter values resulting from fit. First one is 'A', second one is 'B'.
    % Use coeffnames(fit_msd_line) to find out order of parameters.
    fit_msd_line_offset = fit_param_values(1); % offset from fit in pix^2.
    fit_msd_line_slope = fit_param_values(2); % slope from fit in pix^2/s.
    fit_msd_line_rsq = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_msd_line,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    fit_msd_line_offset_stDev = errorSTDEV(1);
    fit_msd_line_slope_stDev = errorSTDEV(2);
catch ME1 % catch error. If there was an error:
    fit_msd_line = [0]; % needed for plot later.
    fit_msd_line_offset = []; 
    fit_msd_line_offset_stDev = []; 
    fit_msd_line_slope = []; 
    fit_msd_line_slope_stDev = []; 
    fit_msd_line_rsq = []; 
end

disp(' ') % empty line
disp('Fit for msd versus delta-time to a line:') 
disp([' fit_msd_line_offset = ',num2str(fit_msd_line_offset),' +- ',num2str(fit_msd_line_offset_stDev),';   fit_msd_line_slope = ',num2str(fit_msd_line_slope),' +- ',num2str(fit_msd_line_slope_stDev)]) 
disp([' fit_msd_line_rsq = ',num2str(fit_msd_line_rsq)]) % r squared of fit.

% Mobility results for msd vs delta-time (numbers):
results_mobility.lengthMsdVector = length(ydata_msd); % no. of points in msd vector with error < 150%.
results_mobility.fit_msd_line_offset = fit_msd_line_offset; 
results_mobility.fit_msd_line_offset_stDev = fit_msd_line_offset_stDev;
results_mobility.fit_msd_line_slope = fit_msd_line_slope;
results_mobility.fit_msd_line_slope_stDev = fit_msd_line_slope_stDev;
results_mobility.fit_msd_line_rsq = fit_msd_line_rsq;

% -----------------------
% Fit the MSD versus delta-time to a SATURATING CURVE (confined mobility):
% As before, only fit data points with relative error < 150%, i.e., xdata_msd and ydata_msd:
% Note: min. no. of points in track must be at least 6 to get at least 4
% points in msd so that a fit with 3 params works ok, otherwise, it fails.
saturating_function = fittype('C+A*(1-exp(-t/B))','independent','t'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% CHECK: values of guesses:
guess_A = 3; % saturation limit guess (pix^2).
guess_B = 0.2; % time-constant guess in seconds.
guess_C = 0; % offset guess (pix^2).
options.StartPoint = [guess_A guess_B guess_C];
% Error control: catch error message to avoid exiting the whole function:
try
    [fit_msd_conf gof] = fit(xdata_msd,ydata_msd,saturating_function,options);
    fit_param_values = coeffvalues(fit_msd_conf); % parameter values resulting from fit. First one is 'A', second one is 'B'.
    % Use coeffnames(fit_msd_conf) to find out order of parameters.
    fit_msd_conf_limit = fit_param_values(1); 
    fit_msd_conf_timeconst = fit_param_values(2); 
    fit_msd_conf_offset = fit_param_values(3);
    fit_msd_conf_rsq = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_msd_conf,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    fit_msd_conf_limit_stDev = errorSTDEV(1);
    fit_msd_conf_timeconst_stDev = errorSTDEV(2);
    fit_msd_conf_offset_stDev = errorSTDEV(3);
catch ME1 % MessageError1. % If fit fails:
    try
        % CHECK: second guesses:
        options.StartPoint = [max(ydata_msd) guess_B guess_C]; % use a different guess for the saturating limit.
        [fit_msd_conf gof] = fit(xdata_msd,ydata_msd,saturating_function,options);
        fit_param_values = coeffvalues(fit_msd_conf); % parameter values resulting from fit. First one is 'A', second one is 'B'.
        fit_msd_conf_limit = fit_param_values(1); % offset from fit. This is the initial intensity, I at t_origin or tstart.
        fit_msd_conf_timeconst = fit_param_values(2); % slope from fit.
        fit_msd_conf_offset = fit_param_values(3);
        fit_msd_conf_rsq = gof.rsquare; % rsquare coefficient of fit.
        errors = confint(fit_msd_conf,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
        errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
        fit_msd_conf_limit_stDev = errorSTDEV(1);
        fit_msd_conf_timeconst_stDev = errorSTDEV(2);
        fit_msd_conf_offset_stDev = errorSTDEV(3);
    catch ME2 % if second try of fit fails too, set all results to zero:
        disp('Fit of msd versus Delta-time to a saturating curve failed.');
        
        fit_msd_conf_limit = 0;
        fit_msd_conf_timeconst = 0;
        fit_msd_conf_offset = 0;
        fit_msd_conf_rsq = 0; % this is safe for later mobility flags too.
        fit_msd_conf_limit_stDev = 0;
        fit_msd_conf_timeconst_stDev = 0;
        fit_msd_conf_offset_stDev = 0;
        
        disp(ME2.message);
        ME2 = addCause(ME2, ME1);
        % rethrow(ME2)
    end
end

disp(' ') % empty line
disp('Fit for msd versus delta-time to a saturating curve:') 
disp([' fit_msd_conf_limit = ',num2str(fit_msd_conf_limit),' +- ',num2str(fit_msd_conf_limit_stDev),';   fit_msd_conf_timeconst = ',num2str(fit_msd_conf_timeconst),' +- ',num2str(fit_msd_conf_timeconst_stDev),';   fit_msd_conf_offset = ',num2str(fit_msd_conf_offset),' +- ',num2str(fit_msd_conf_offset_stDev)]) 
disp([' fit_msd_conf_rsq = ',num2str(fit_msd_conf_rsq)]) % r squared of fit.

% Mobility results for msd vs delta-time (numbers):
results_mobility.fit_msd_conf_limit = fit_msd_conf_limit; 
results_mobility.fit_msd_conf_limit_stDev = fit_msd_conf_limit_stDev;
results_mobility.fit_msd_conf_timeconst = fit_msd_conf_timeconst;
results_mobility.fit_msd_conf_timeconst_stDev = fit_msd_conf_timeconst_stDev;
results_mobility.fit_msd_conf_offset = fit_msd_conf_offset;
results_mobility.fit_msd_conf_offset_stDev = fit_msd_conf_offset_stDev;
results_mobility.fit_msd_conf_rsq = fit_msd_conf_rsq;

% ----------------
% Plot only msd data which has relative error < 150% :
subplot(3,3,5);
errorbar(xdata_msd,ydata_msd,halfErrorBars,halfErrorBars,'.b'); 
xlim([0 max(xdata_msd)]); 
ylim([-0.5*max(ydata_msd) 2.5*max(ydata_msd)]);  
hold on;
% Plot msd fits:
plot(fit_msd_line,'k'); % plot linear fit to msd data as a black line.
try % error control, if fit to saturating curve failed, do not plot:s
    plot(fit_msd_conf,'r'); % plot saturating-curve fit to msd data as a red line.
catch ME3 % catch error to avoid getting out of the entire function which is runnings
end
legend('hide');
xlabel('\Deltat (s)'); 
ylabel('msd (pix^2)'); 
hold off;

% multifigure continues to be filled in later...

% -------------
% Get values of size of 1D confinement region and Diffusion coefficient in
% meaningful units.
% For the "Brownian" trajectories: from the slope of msd vs Delta-t I get the
% 1D Diffusion coefficient as slope/2 and use the pixel size to convert to
% microns^2/s. 
% The slope from the fit of msd vs deltat-time is in pix^2/s.
diffusion1D_coeff_micronsSqrdPerSec = (fit_msd_line_slope/2)*(pixelsize_nm/1000)^2; % 1D diffusion coeff in microns^2/s.
diffusion1D_coeff_error = (fit_msd_line_slope_stDev/2)*(pixelsize_nm/1000)^2; % error or standard deviation of the 1D diffusion coeff in microns^2/s.
disp(' ') % empty line
disp(['Diffusion1D_coeff_micronsSqrdPerSec = ',num2str(diffusion1D_coeff_micronsSqrdPerSec),' +- ',num2str(diffusion1D_coeff_error)])

% For the "confined" trajectories, I assume confined diffusion in 1D to a
% segment of length L. The msd saturates to a value of (L^2)/6. I use the
% pixel size to convert to nm. From the fit, the value the msd saturates to
% is (A) fit_msd_conf_limit (pixels^2). 
size_1Dconfin_nm = pixelsize_nm*sqrt(6*fit_msd_conf_limit); % size of 1D confinement region in nm (L).
% L = pixtonm*sqrt(6)*sqrt(A), where A =fit_msd_conf_limit and L = size_1Dconfin_nm. 
% Propagating errors we have that error_L = (L/2)*(error_A/A):
size_1Dconfin_nm_error = (size_1Dconfin_nm/2)*(fit_msd_conf_limit_stDev/fit_msd_conf_limit); % size of 1D confinement region in nm.
disp(' ') % empty line
disp(['size_1Dconfin_nm = ',num2str(size_1Dconfin_nm),' +- ',num2str(size_1Dconfin_nm_error)])

% ----------------------
% Mobility flags:
if fit_msd_line_rsq > mobility_rsq_limit
    Brownian_flag = 1;
else 
    Brownian_flag = 0;
end

if (fit_msd_conf_rsq > mobility_rsq_limit) && (fit_msd_conf_timeconst < max_conf_timeconst) && (fit_msd_conf_timeconst < 1.5*traj_duration)
    % Reject fits to confined trajectory if timeconstant of fit is >= max_conf_timeconst
    % seconds (that is really a linear fit):
    Confined_flag = 1;
else
    Confined_flag = 0;
end

% (Brownian_flag == 1 && Confined_flag == 1) ||
if  Brownian_flag == 0 && Confined_flag == 0
    OtherMobility_flag = 1; % this results in mobility flags [1 1 1] or [0 0 1]
else
    OtherMobility_flag = 0;
end

mobility_flags = [Brownian_flag Confined_flag OtherMobility_flag];

disp(' '); % empty line.
disp('Mobility flags [Brownian_flag  Confined_flag  OtherMobility_flag] are now: ');
disp(mobility_flags);
% CHECK: give user option to override Mobility Flags or not:
% Give user a chance to override mobility flags:
% agree_mobility_flags = input('Do you agree with previous mobility_flags? (1 for yes, 0 for no): '); % request user input
agree_mobility_flags = 1;
if agree_mobility_flags == 0 % if user does not agree with mobility flags:
    mobility_flags = input('Enter new mobility_flags as [traj_Brownian  traj_confined  traj_other] (ones and zeros):'); % request user input.
    Brownian_flag = mobility_flags(1);
    Confined_flag = mobility_flags(2);
    OtherMobility_flag = mobility_flags(3);
end

% ------------------------
% Add to mobility results (numbers):
results_mobility.diffusion1D_coeff_micronsSqrdPerSec = diffusion1D_coeff_micronsSqrdPerSec;
results_mobility.diffusion1D_coeff_error = diffusion1D_coeff_error;
results_mobility.size_1Dconfin_nm = size_1Dconfin_nm;
results_mobility.size_1Dconfin_nm_error = size_1Dconfin_nm_error;
results_mobility.Brownian_flag = Brownian_flag;
results_mobility.Confined_flag = Confined_flag;
results_mobility.OtherMobility_flag = OtherMobility_flag;



%% Calculate Intensity pair-wise differences and their probability distribution:

if use_filtered == 1  % use filtered data.
    I_to_use = filtered_I; 
    % otherwise, I_to_use stays as it was, use unfiltered data. 
end

% Calculate distribution of pair-wise differences (see fullPwD.m):
% nbins: see PARAMETERS section: number of bins to use.
[binCentres freqCounts] = fullPwD(I_to_use,nbins,0); % Calculate distribution of pair-wise differences.

% Select only positive and only negative values of pairwise intensity differences (horiz axis
% in the previous histogram) and plot them as a solid line:
pos_positives = find(binCentres > 0); % positions of positive bin-centre values.
pos_negatives = find(binCentres < 0); % positions of negative bin-centre values.

% Error control:
if isempty(pos_negatives) % this is a fast dodgy way of sorting this out for the moment.
    % pos_negatives = 1;
    disp('!!!  No negative side of histogram !!!');
    binCentres_neg = [0];
    freqCounts_neg = [0];
else
    binCentres_neg = binCentres(pos_negatives);
    freqCounts_neg = freqCounts(pos_negatives);
end

if isempty(pos_positives) % this is a fast dodgy way of sorting this out for the moment.
    % pos_positives = 1;
    disp('!!! No positive side of histogram !!!');
    binCentres_pos = [0];
    freqCounts_pos = [0];
else
    binCentres_pos = binCentres(pos_positives);
    freqCounts_pos = freqCounts(pos_positives);
end

% Plot distributions of pair-wise differences:
% figure(2)
subplot(3,3,3);
bar(binCentres,freqCounts,'c'); % plot a bar graph of the full histogram.
hold on; 
plot(binCentres_neg,freqCounts_neg,'k-','LineWidth',1); % superimpose continuous line only for positive side.
legend('All intensity jumps','Intensity drops','Location','Best');
xlabel('Intensity pair-wise differences'); 
ylabel('frequency');
ylim([0 1.1*max(freqCounts)]); % re-scale vertical axis.
hold off;


%% Calculate spectrum (Fourier transform) of pair-wise difference
% distribution (only "Negative" values) and find peaks in it:
% "Negative" means I(i+1)<I(i), i.e., drops in intensity as time increases.

% Calculate and plot power spectrum (ps) of distribution of pair-wise
% differences and peaks found in it, in order of increasing spectral power
% (see FourierAndFindPeaks.m).
% Use only the negative binCentres values, binCentres_neg, to compute the
% spectrum. 

[ps_x_neg ps_y_neg ps_peaks_x_neg ps_peaks_y_neg] = FourierAndFindPeaks(binCentres_neg,freqCounts_neg,0,0);
% last input = 0 means figure is not displayed within function "FourierAndFindPeaks.m".

% Plot spectral power versus intensity_step values.
% Plot single-sided power spectrum (only of positive intensity pair-wise differences):
% figure(3)
subplot(3,3,6);
% Power spectrum of prob distrib. of I-pair-wise differences:
    plot(ps_x_neg,ps_y_neg,'b-','LineWidth',1.5);
    xlim([0 x_limit_spectr]);
    title({['Traj ',num2str(n_traj),' - PowerSpectr (prob distrib I drops)']});
    xlabel('Intensity period step')
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow = min(length(ps_peaks_x_neg),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow
        label_to_print = strcat('\leftarrow',num2str(ps_peaks_x_neg(i),3)); % print 3 significant figures.
        text(ps_peaks_x_neg(i),ps_peaks_y_neg(i),label_to_print,'FontSize',10);
    end


%% Calculate spectrum (Fourier transform) of POSITIVE side of histogram of pair-wise differences:
% and find peaks in it. Plot. "Positive" means I(i+1)>I(i), i.e., increasing intensity jumps.

% Now use only the positive binCentres values, binCentres_neg, and look at
% their spectrum:
[ps_x_pos ps_y_pos ps_peaks_x_pos ps_peaks_y_pos] = FourierAndFindPeaks(binCentres_pos,freqCounts_pos,0,0);
% last input = 0 means figure is not displayed within function "FourierAndFindPeaks.m".

% Plot spectral power versus intensity_step values (only of positive
% intensity pair-wise differences):
subplot(3,3,9);
% Power spectrum of prob distrib. of I-pair-wise differences:
    plot(ps_x_pos,ps_y_pos,'b-','LineWidth',1.5);
    xlim([0 x_limit_spectr]);
    title({['Traj ',num2str(n_traj),' - PowerSpectr (prob distrib of I increases)']});
    xlabel('Intensity period step') 
    ylabel('|Fourier Transform|^2')
    % Mark highest found peaks on previous graph:
    numOfPeaksToShow_pos = min(length(ps_peaks_x_pos),8); % show only the highest 8 peaks or less if there are less than 8 peaks found.
    for i=1:numOfPeaksToShow_pos
        label_to_print_pos = strcat('\leftarrow',num2str(ps_peaks_x_pos(i),3)); % print 3 significant figures.
        text(ps_peaks_x_pos(i),ps_peaks_y_pos(i),label_to_print_pos,'FontSize',10);
    end



%% Plot frame average with all trajectory points:

% Display frame average and all trajectory positions overlaid on it:
subplot(3,3,8);
imshow(frame_avg,[],'Border','tight'); % display frame average scaled between its min and max values ([]).
hold on;
% x_values = analysedAllTraj(n_traj).xvalues0; % list of original x centres of spots in trajectory n_traj.
% y_values = analysedAllTraj(n_traj).yvalues0; % list of original y centres of spots in trajectory n_traj.

if mean(y_values) < (frame_Ysize/2)
    circleColour_string = 'r'; % plot circles in red if in top half of image.
else
    circleColour_string = 'g'; % plot circles in red if in bottom half of image.
end

for k = 1:length(x_values)
    plot(x_values(k),y_values(k),'.-','Color',circleColour_string,'MarkerSize',5) % plot track positions in green.
end
hold off;



%% Continue adding results to figure:

% % Reminder:
% firstFrameInTraj = analysedAllTraj(n_traj).frame(1); % first frame in trajectory.
% lastFrameInTraj = analysedAllTraj(n_traj).frame(end); % last frame in trajectory.
% traj_duration = tsamp*(lastFrameInTraj-firstFrameInTraj); % duration of track in seconds.
% nPointsInTrack = analysedAllTraj(n_traj).numel;

subplot(3,3,7); axis off; 
% text(-0.07,1,trajXlsPath,'Interpreter','none'); % display path of image sequence at position (0,1) of the axes (when there is no specific plot, axes go from 0 to 1).

str1(1) = {['Plots shown are for TRAJECTORY number ',num2str(n_traj)]};
str1(2) = {['No. of data points in this trajectory: ',num2str(nPointsInTrack)]};
str1(3) = {['First frame in this trajectory is:  ',num2str(firstFrameInTraj)]};
str1(4) = {['Time origin for sequence is frame no. ',num2str(start_frame)]};
str1(5) = {['Time between frames in seconds is: ',num2str(tsamp)]};
str1(6) = {['No. trajectories analysed: ',num2str(length(analysedAllTraj))]};
str1(7) = {['original TRAJECTORY no.: ',num2str(analysedAllTraj(n_traj).trajNumber)]};

% str1(8) = {['Fit I vs time to Ioffset+I0*exp(-(t-tstart)/tau), with tstart=const= ',num2str(t_origin),'.']};
% str1(9) = {['I0 = ',num2str(I0_fit),' \pm ',num2str(stDev_I0)]};
% str1(10) = {['tau = ',num2str(tau_fit),' \pm ',num2str(stDev_tau),' s.']};
% str1(11) = {['Ioffset = ',num2str(Ioffset_fit),' \pm ',num2str(stDev_Ioffset)]};
% str1(12) = {['r square of fit = ',num2str(rsq_fit_I)]};
% str1(13) = {['Window size for C-K filter (points) = ',num2str(Wfilter)]};
% str1(14) = {['Weighting exponent for C-K filter = ',num2str(Rfilter)]};
% str1(15) = {['Use filtered intensity trace? ',num2str(use_filtered)]};

% Print intensity step sizes (found peaks) only for positive steps:

% % Now print only the highest peaks:
% if ~isempty(ps_peaks_x_neg)  % if peaks were found
% str1(16) = {['Found I steps = ',num2str(ps_peaks_x_neg',4)]};
% str1(17) = {['Estimated no. molecules = ',num2str(ps_peaks_y_neg',4)]};
% end

text(0,0.5,str1)

% % Display info for Intensity vs time fit and for Chung-Kennedy filter parameters: 
% subplot(3,3,8); axis off;
% str2(1) = {['Fit I vs time to Ioffset+I0*exp(-(t-tstart)/tau), with tstart=const= ',num2str(t_origin),'.']};
% str2(2) = {['I0 = ',num2str(I0_fit),' \pm ',num2str(stDev_I0)]};
% str2(3) = {['tau = ',num2str(tau_fit),' \pm ',num2str(stDev_tau),' s.']};
% str2(4) = {['Ioffset = ',num2str(Ioffset_fit),' \pm ',num2str(stDev_Ioffset)]};
% str2(5) = {['r square of fit = ',num2str(rsq_fit_I)]};
% str2(6) = {['Window size for C-K filter (points) = ',num2str(Wfilter)]};
% str2(7) = {['Weighting exponent for C-K filter = ',num2str(Rfilter)]};
% str2(8) = {['Use filtered intensity trace? ',num2str(use_filtered)]};
% text(0,0.5,str2)

% % str3(1) = {'FIT RESULTS for msd vs \Deltat:'};
% % str3(2) = {'Fit to equation:  '};
% % str3(3) = {' '};
% str3(1) = {['Found I steps = ',num2str(Istep_peaks)]};
% str3(2) = {['Estimated no. molecules = ',num2str(num_molecs)]};
% % str3(6) = {['r square of fit = ',num2str(rsq_fit_I)]};
% text(0,0.5,str3)
% Note: I can write in LateX form in the graph labels... t_{abs} appears on
% graph with 'abs' as a subindex and \Delta gives me the greek Delta
% character...!!!


%% ADD FLAGS:
% Add flags to first element of cell array of results:

n_trackInfo = 1; % cell element in which to save track info.
figure(h1); % bring figure to front.

% Flag indicating if all points in track are equally spaced or if there are
% jumps in track (due to linking segments up to a certain no. of frames away):
if traj_duration > tsamp*(nPointsInTrack-1)
    track_with_jumps_flag = 1;
else
    track_with_jumps_flag = 0;
end
processedTraj{n_trackInfo}.track_with_jumps_flag = track_with_jumps_flag; 

% Good-tracking flag:
% only files which show "good tracking" (on video) in showManyTrajAnalysis.m are analysed, and their results saved.
good_tracking_flag = 1;

% Dividing-cell flag:
% disp(' ');
% dividingCell_flag = input('Is this a dividing cell? (1 for yes, 0 for no): '); % request user input
dividingCell_flag = 0;

% Flag trajectories as "short" (see PARAMETERS section):
if nPointsInTrack < max_NumFramesForShortTrack
    short_track_flag = 1;
else
    short_track_flag = 0;
end

% Flag trajectories as "long" or "very long" (see PARAMETERS section):
if nPointsInTrack < min_NumFramesForLongTrack
    long_track_flag = 0;
    very_long_track_flag = 0;
else
    if nPointsInTrack < min_NumFramesForVeryLongTrack
        long_track_flag = 1;
        very_long_track_flag = 0;
    else
        long_track_flag = 0;
        very_long_track_flag = 1;
    end    
end

% FLAG trajectories which are GOOD FOR getting INITIAL INTENSITY:
% If the first frame in this trajectory is less than
% "max_framesAwayFromTimeOrigin" frames away from the time-origin
% frame for this image-sequence, then flag as good to extrapolate initial
% intensity (from the linear fit to the first points in track):
if (firstFrameInTraj-start_frame) <= max_framesAwayFromTimeOrigin
    track_closeToTimeOrigin_flag = 1;
else
    track_closeToTimeOrigin_flag = 0;
end

% FLAG trajectories which show GOOD EXPONENTIAL FITS of I vs time:
% Flag as good exp fit if timeconstant (tau) from fit no larger than the total
% track duration, if relative error of tau below max_tau_relativeError_forGoodExpFit (%) and if rsq of fit > min_rsq_forGoodExpFit:
% For exponential fit with no offset:
try % error control in case the fits failed:
    if (tau_fit <= traj_duration) && (100*stDev_tau/tau_fit < max_tau_relativeError_forGoodExpFit) && ...
            (rsq_fit_I > min_rsq_forGoodExpFit)
        goodExpFit_flag = 1;
    else
        goodExpFit_flag = 0;
    end
catch
    goodExpFit_flag = [];
end

% For exponential fit with offset:
try
    if (tau_fit_wo <= traj_duration) && (100*stDev_tau_wo/tau_fit_wo < max_tau_relativeError_forGoodExpFit) && ...
            (rsq_fit_I_wo > min_rsq_forGoodExpFit)
        goodExpFit_wo_flag = 1;
    else
        goodExpFit_wo_flag = 0;
    end
catch
    goodExpFit_wo_flag = [];
end
% TOP OR BOTTOM of image (Green (bottom) or Red (top) channels):
% Image size is frame_Xsize by frame_Ysize in pixels, so if
% average y value is below half frame_Ysize then we are at top, otherwise, at bottom of image:
if mean(analysedAllTraj(n_traj).yvalues0) < (frame_Ysize/2)
    TopOrBottom_flag = 'top';
else
    TopOrBottom_flag = 'bottom';
end

% Flag tracks with low enough intensity levels above background (to use
% later to try and get stoichiometry from intensity steps in track).
% Use a different I limit for top and bottom channels (GFP, mCherry):
if strcmp(TopOrBottom_flag,'bottom') % string compare: true or false.
    pos_below_Ilimit = find(IforFit<lowEnoughI_limit_bottom); % find positions in vector.
else % if TopOrBottom_flag is 'top'
    pos_below_Ilimit = find(IforFit<lowEnoughI_limit_top); % find positions in vector.
end
I_belowLimit = IforFit(pos_below_Ilimit); % select only those intensities below Ilimit (lowEnoughI_limit_top or lowEnoughI_limit_bottom).
t_for_IbelowLimit = tforFit(pos_below_Ilimit); % select the accompanying times;
% Add flag for low enough intensities above background level:
if isempty(I_belowLimit)
    lowEnoughI_flag = 0; % No intensities below limit.
else
    lowEnoughI_flag = 1; % Some intensities below limit.
end

% Save flags into TRACK INFO results:
processedTraj{n_trackInfo}.good_tracking_flag = good_tracking_flag; 
processedTraj{n_trackInfo}.TopOrBottom = TopOrBottom_flag;
processedTraj{n_trackInfo}.dividingCell_flag = dividingCell_flag;
processedTraj{n_trackInfo}.short_track_flag = short_track_flag; 
processedTraj{n_trackInfo}.long_track_flag = long_track_flag;
processedTraj{n_trackInfo}.very_long_track_flag = very_long_track_flag;
processedTraj{n_trackInfo}.track_closeToTimeOrigin_flag = track_closeToTimeOrigin_flag;
processedTraj{n_trackInfo}.lowEnoughI_flag = lowEnoughI_flag;
processedTraj{n_trackInfo}.goodExpFit_flag = goodExpFit_flag;
processedTraj{n_trackInfo}.goodExpFit_wo_flag = goodExpFit_wo_flag;

% Display flags on screen:
disp('Trajectory flags are now: ');
disp(processedTraj{n_trackInfo}); 
figure(h1); % bring figure to front.
% numbers, first column is fieldnames of element {n_trackInfo}, second one is values.
% CHECK: give user chance to override TRAJECTORY FLAGS or not:
% Give user a chance to override trajectory flags:
% agree_track_flags = input('Do you agree with previous Trajectory Flags? (1 for yes, 0 for no): '); % request user input
agree_track_flags = 1;
if agree_track_flags ~= 1 % if user does not agree with track flags:
    track_flags = input('Enter new Track Flags as [good_tracking  track_closeToTimeOrigin  goodExpFit  goodExpFit_withOffset] (1 or 0 each):'); % request user input.
    good_tracking_flag = track_flags(1);
    track_closeToTimeOrigin_flag = track_flags(4);
    goodExpFit_flag = track_flags(5);
    goodExpFit_wo_flag = track_flags(6);
    
    % Save new flags into TRACK INFO results:
    processedTraj{n_trackInfo}.good_tracking_flag = good_tracking_flag;
    processedTraj{n_trackInfo}.track_closeToTimeOrigin_flag = track_closeToTimeOrigin_flag;
    processedTraj{n_trackInfo}.goodExpFit_flag = goodExpFit_flag;
    processedTraj{n_trackInfo}.goodExpFit_wo_flag = goodExpFit_wo_flag;
end

% Save deciding criteria for flags into TRACK INFO results:
processedTraj{n_trackInfo}.max_NumFramesForShortTrack = max_NumFramesForShortTrack; 
processedTraj{n_trackInfo}.min_NumFramesForLongTrack = min_NumFramesForLongTrack; 
processedTraj{n_trackInfo}.min_NumFramesForVeryLongTrack = min_NumFramesForVeryLongTrack;
processedTraj{n_trackInfo}.max_framesAwayFromTimeOrigin = max_framesAwayFromTimeOrigin;
processedTraj{n_trackInfo}.lowEnoughI_limit_top = lowEnoughI_limit_top;
processedTraj{n_trackInfo}.lowEnoughI_limit_bottom = lowEnoughI_limit_bottom;
processedTraj{n_trackInfo}.min_rsq_forGoodExpFit = min_rsq_forGoodExpFit;
processedTraj{n_trackInfo}.max_tau_relativeError_forGoodExpFit = max_tau_relativeError_forGoodExpFit;


%% Save graphical analysis results (last figure):

% figName = strcat('figResultsTraj',num2str(n_traj)); % name of figure file to save to.
% set(gcf, 'PaperPositionMode', 'auto')  % Use screen size.
% print('-dpng',figName)  % add -r200 after -dpng to control resolution.

% Move into folder previously created to save traj analysis results and
% SAVE current FIGURE as a .png (at screen size) and then close the figure window:
% Result figure files are saved in a new folder within the same directory where the input excel file was.
figName = strcat(new_folder_name,'_figResults_traj',num2str(n_traj)); % name of figure file to save to. 
saveFigurePNG(new_folder_name,figName); % See saveFigurePNG.m.


%% Output, analysis results:
% Prepare result data to save them later.

% TRACK INFO continued (numbers):
% First element in output cell array contains general info about the image and the trajectory.
% Output result: 
% n_trackInfo = 1; % cell element in which to save track info. See above.
processedTraj{n_trackInfo}.XLSpath = trajXlsPath;
processedTraj{n_trackInfo}.minNumPointsInTraj = analysedAllTraj(n_traj).minNumPointsInTraj; % min no. of points in trajectory for it to be analysed.
processedTraj{n_trackInfo}.ImageSizeHorizPix = frame_Xsize;
processedTraj{n_trackInfo}.ImageSizeVerticPix = frame_Ysize;
processedTraj{n_trackInfo}.TrajNum = n_traj; % numbers of trajectories analysed (long enough) as chosen in analyseTraj.m.
processedTraj{n_trackInfo}.OriginalTrajNum = analysedAllTraj(n_traj).trajNumber; % Original traj number from all those tracked initially.
processedTraj{n_trackInfo}.FirstTrajFrame = firstFrameInTraj; % First frame in this trajectory.
processedTraj{n_trackInfo}.LastTrajFrame = lastFrameInTraj; % Last frame in this trajectory.
processedTraj{n_trackInfo}.NumDataPoints = analysedAllTraj(n_traj).numel; % no. of points in track.
processedTraj{n_trackInfo}.TrajStartTime = tsamp*analysedAllTraj(n_traj).frame(1); % absolute start time of trajectory (s).
processedTraj{n_trackInfo}.TrajEndTime = tsamp*analysedAllTraj(n_traj).frame(end); % absolute end time of trajectory (s).
processedTraj{n_trackInfo}.TrajDuration = traj_duration; % duration of track in seconds.
processedTraj{n_trackInfo}.TimeBetweenFrames = tsamp;
processedTraj{n_trackInfo}.FrameForTimeOrigin = start_frame; % time origin (as frame number) for this whole image sequence.
processedTraj{n_trackInfo}.AbsTimeOrigin = t_origin;
processedTraj{n_trackInfo}.pixelsize_nm = pixelsize_nm; % pixel size in nm.
processedTraj{n_trackInfo}.Track_meanX_0 = mean(x_values); % mean x position (original) of track on image.
processedTraj{n_trackInfo}.Track_meanY_0 = mean(y_values); % mean y position (original) of track on image.
processedTraj{n_trackInfo}.Track_meanX_cellOrigin = mean(x_values_2); % mean x-x_com position of track with respect to cell centre of mass.
processedTraj{n_trackInfo}.Track_meanY_cellOrigin = mean(y_values_2); % mean y-y_com position of track with respect to cell centre of mass.
processedTraj{n_trackInfo}.Track_meanXvalue_cellCoords = mean(x_values_cell); % mean x' position of track in local cell coordinates.
processedTraj{n_trackInfo}.Track_meanYvalue_cellCoords = mean(y_values_cell); % mean y' position of track in local cell coordinates.
% Parameters for Chung-Kennedy "CK" filter:
processedTraj{n_trackInfo}.WindowWidthCKfilter = Wfilter;
processedTraj{n_trackInfo}.WeightExponentCKfilter = Rfilter;
processedTraj{n_trackInfo}.useFiltered = use_filtered; % use filtered intensity trace or not.
% Parameter for number of bins in histogram of intensity pair-wise differences:
processedTraj{n_trackInfo}.numBins = nbins; 
processedTraj{n_trackInfo}.meanSpotWidth = meanSpotWidth; % Average fitted Gaussian width of tracked spot, in pixels.

% -----------------
% TRACK DATA (vectors):
% Next element of results cell array contains track data:
n_trackData = 2; % cell element in which to save track info.
processedTraj{n_trackData}.frame = analysedAllTraj(n_traj).frame; % vectors.
processedTraj{n_trackData}.timeabs = analysedAllTraj(n_traj).timeabs;
processedTraj{n_trackData}.intensity = analysedAllTraj(n_traj).intensity; % This is IspTot, integrated spot intensity.
processedTraj{n_trackData}.filtered_I = filtered_I; % Chung-Kennedy Filtered trace.
processedTraj{n_trackData}.I0_IspotFit = analysedAllTraj(n_traj).I0_IspotFit;
processedTraj{n_trackData}.sigma_IspotFit = analysedAllTraj(n_traj).sigma_IspotFit;
processedTraj{n_trackData}.xvalues0 = analysedAllTraj(n_traj).xvalues0; % original position on image.
processedTraj{n_trackData}.yvalues0 = analysedAllTraj(n_traj).yvalues0; % original position on image.
processedTraj{n_trackData}.xvalues = analysedAllTraj(n_traj).xvalues; % position with respect to 1st point in track.
processedTraj{n_trackData}.yvalues = analysedAllTraj(n_traj).yvalues; % position with respect to 1st point in track.
% trajectory data in the cell-coordinate system:
processedTraj{n_trackData}.x_values_cell = x_values_cell; % vector, traj in local cell coordinates.
processedTraj{n_trackData}.y_values_cell = y_values_cell; % vector, traj in local cell coordinates.
% Others:
processedTraj{n_trackData}.SNR = analysedAllTraj(n_traj).SNR; % Signal to noise ratio.
processedTraj{n_trackData}.bg_noise_offset_afterBGsubtract = analysedAllTraj(n_traj).bg_noise_offset_afterBGsubtract; % Offset background noise per pixel on average after background subtraction.
processedTraj{n_trackData}.BgNoiseStd = analysedAllTraj(n_traj).BgNoiseStd; % Background noise standard deviation.
processedTraj{n_trackData}.IbgAvg = analysedAllTraj(n_traj).IbgAvg; % Intensity in background region, on average and per pixel.
processedTraj{n_trackData}.IinnerTot = analysedAllTraj(n_traj).IinnerTot; % Total intensity in circular signal mask (no bgnd subtraction).
processedTraj{n_trackData}.rsqFit  = analysedAllTraj(n_traj).rsqFit; % r-square of Gaussian fit of spot intensity vs distance to spot centre for all pixels in spot signal circular mask.

% -----------------
% RESULTS I FITS (numbers):
% Next element of results cell array (function output) contains results 
% from exponential fits to intensity vs time (with and without offset) and
% from fit of first few intensity points to a line to get initial
% intensity:
n_results_Ifits = 3; % cell element in which to save track info.
processedTraj{n_results_Ifits} = results_I_fits;


% -----------------
% MSD RESULTS (vectors):
% Next element of results cell array contains MSD results:
n_MSD = 4;
% Msd data points for which error < 150% :
processedTraj{n_MSD}.deltaTime = xdata_msd; % vectors.
processedTraj{n_MSD}.msd = ydata_msd;
processedTraj{n_MSD}.errorMsd = halfErrorBars;
processedTraj{n_MSD}.errorMsdRelPercent = msd_relative_errors;

% % The original msd data (all points):
% processedTraj{n_MSD}.deltaTime = analysedAllTraj(n_traj).deltaTime; % vectors.
% processedTraj{n_MSD}.msd = analysedAllTraj(n_traj).msd;
% processedTraj{n_MSD}.errorMsd = analysedAllTraj(n_traj).errorMsd;
% processedTraj{n_MSD}.errorMsdRelPercent = analysedAllTraj(n_traj).errorMsdRelPercent;


% -----------------
% MOBILITY RESULTS from fits (numbers):
n_mobility = 5;
processedTraj{n_mobility} = results_mobility;

% -----------------
% HISTOGRAMS (vectors):
% Histograms of Intensity pairwise differences: 
n_hist = 6;
processedTraj{n_hist}.BinCentresPos = binCentres_pos';
processedTraj{n_hist}.FreqCountsPos = freqCounts_pos';

% -----------------
% FOURIER (vectors): 
% Spectrum (Fourier transform) of the previous histogram:
n_Fourier = 7;
processedTraj{n_Fourier}.IstepAxis = ps_x_neg;
processedTraj{n_Fourier}.PowerSpectrum = ps_y_neg;

% -----------------
% FOURIER PEAKS (vectors):
% Periodicity of intensity steps from photobleaching steps of intensity (vectors):
n_FourierPeaks = 8;
if ~isempty(ps_peaks_x_neg)  % if peaks are found (at least one)
processedTraj{n_FourierPeaks}.Isteps = ps_peaks_x_neg; % vector with intensity steps, i.e., found peaks in spectrum.
processedTraj{n_FourierPeaks}.power_Fourier = ps_peaks_y_neg; % vector with Fourier spectral powers, i.e., height of found peaks in spectrum.
else
processedTraj{n_FourierPeaks}.Isteps = []; 
processedTraj{n_FourierPeaks}.power_Fourier = []; 
end

% -----------------
% LOCAL CELL COORDINATES (numbers):
% Results for local cell coordinates:
n_localCellCoords = 9;
processedTraj{n_localCellCoords} = results_cell_coords;
% note that the angle in the image plane has oposite sign to the angle in
% an x-y plot, because image plane has reflected y coord which grows downwards.
% -----------------------------------------

% -----------------
% I BELOW LIMIT (vectors): 
% Save to another sheet only the intensity values below the limits lowEnoughI_limit_top or lowEnoughI_limit_bottom:
n_IbelowLimit = 10;
processedTraj{n_IbelowLimit}.t_for_IbelowLimit = t_for_IbelowLimit;
processedTraj{n_IbelowLimit}.I_belowLimit = I_belowLimit;


%% Save excel file with trajectory results:

% Move into folder previously created to save traj analysis results:
cd(new_folder_name); % move into that directory.

% Save results in another excel file in the same folder where the
% analysis plot .png file has been saved, new_folder_name:
output_filename = strcat(new_folder_name,'_traj',num2str(n_traj),'.xls'); % name of excel file to save to.
dataForSheet1 = [fieldnames(processedTraj{n_trackInfo}) struct2cell(processedTraj{n_trackInfo})]; % numbers, first column is fieldnames of element {n_trackInfo}, second one is values.
dataForSheet2 = [fieldnames(processedTraj{n_trackData})'; num2cell(cell2mat(struct2cell(processedTraj{n_trackData})'))]; % vectors
dataForSheet3 = [fieldnames(processedTraj{n_results_Ifits}) struct2cell(processedTraj{n_results_Ifits})]; % numbers
dataForSheet4 = [fieldnames(processedTraj{n_MSD})'; num2cell(cell2mat(struct2cell(processedTraj{n_MSD})'))]; % vectors
dataForSheet5 = [fieldnames(processedTraj{n_mobility}) struct2cell(processedTraj{n_mobility})]; % numbers
dataForSheet6 = [fieldnames(processedTraj{n_hist})'; num2cell(cell2mat(struct2cell(processedTraj{n_hist})'))]; % vectors
dataForSheet7 = [fieldnames(processedTraj{n_Fourier})'; num2cell(cell2mat(struct2cell(processedTraj{n_Fourier})'))]; % vectors
dataForSheet8 = [fieldnames(processedTraj{n_FourierPeaks})'; num2cell(cell2mat(struct2cell(processedTraj{n_FourierPeaks})'))]; % vectors
dataForSheet9 = [fieldnames(processedTraj{n_localCellCoords}) struct2cell(processedTraj{n_localCellCoords})]; % numbers
dataForSheet10 = [fieldnames(processedTraj{n_IbelowLimit})'; num2cell(cell2mat(struct2cell(processedTraj{n_IbelowLimit})'))]; % vectors

warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.

xlswrite(output_filename,dataForSheet1,'Track info'); % write data to sheet 'Track info' in excel file.
xlswrite(output_filename,dataForSheet2,'Track data'); % write data to sheet 'Track data' in excel file.
xlswrite(output_filename,dataForSheet3,'Results I fits'); % write data to sheet 'Results I fits' in excel file.
xlswrite(output_filename,dataForSheet4,'MSD results'); % write data to sheet 'MSD results' in excel file.
xlswrite(output_filename,dataForSheet5,'Mobility results'); % write data to sheet 'Mobility results' in excel file.
xlswrite(output_filename,dataForSheet6,'Histograms pwd'); % write data to sheet 'Histograms pwd' in excel file.
xlswrite(output_filename,dataForSheet7,'Fourier'); % write data to sheet 'Fourier' in excel file.
xlswrite(output_filename,dataForSheet8,'Fourier Peaks'); % write data to sheet 'Fourier Peaks' in excel file.
xlswrite(output_filename,dataForSheet9,'cell coordinates'); % write data to sheet 'cell coordinates' in excel file.
xlswrite(output_filename,dataForSheet10,'I below Limit'); % write data to sheet 'I below Limit' in excel file.

cd('..'); % go back to previous directory.

% CHECK: skip this pause or not.
% Give time to user to prepare before looking at next track:
% input('Press Enter when ready to look at next track: '); 








