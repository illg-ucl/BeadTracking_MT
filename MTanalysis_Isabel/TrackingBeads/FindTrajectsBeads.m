function bead_results = FindTrajectsBeads(image_label,start_frame,end_frame)
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
% Function to find all bead trajectories in an input image sequence.
%
% INPUTS:
% - image_label: string such as '513', '490', etc., that appears in the
% file name of the image sequence to be analysed, which is in the current folder.
% - start_frame: first frame of the sequence to be analysed.
% - end_frame: last frame of the sequence to be analysed.
% 
% Example of how to call this function: 
% for frames 1 to 10 of image "Heiko_Thu Jun 24 2010_554.sif" in current folder:
% beadResults1 = FindTrajectsBeads('554',1,10); 
% -------------------------------------------------------
% NOTE: for reading .sif files you need 'read_sif_data_direct.m',
% 'GetAndorSifSize', etc. which are currently in path:
% '(...)\myMatlabFiles\openImageSequences\IO_Input\'.
% NOTE: before running this function you should move into the directory which
% contains the image sequence data (labelled by 'image_label').
%
% Reads image sequence data and then for each frame:
% finds candidate bead positions, eliminates coincidences in those candidates (closer than 1-few pixels)
% within all candidates, substracts background,
% finds actual bead centres through cross-correlation method,
% and accepts only bead centres not too close from edge.
% Resulting final accepted bead centres are saved in the structure array bead_results.
% Then link bead positions into trajectory segments: For second frame do differently
% and check only beads in 1st frame and compare to beads in second frame.
% For rest of frames (for frame k): A) first check loose beads (TrajNumber=0) two frames
% ago (k-2) and compare to beadss in current frame (k).
% B) then check all beads in previous frame (k-1) and compare to beads in
% current frame (k).
% To decide on the best asignment of pairs of beads (link), we build up
% matrices of pair-wise distances, ratio of intensities and ratio of sigmas
% of all pairs of beads in the two frames being compared, and take the
% winning asignment as that with the smallest pairwise distance.
% 
% start_frame and end_frames are the frames through which the loop runs to
% find centres of beads.
%
% OUTPUT:
% The output, bead_results is a cell array with two elements:
% bead_results = {params, bead_final};
% The first element in the cell array, bead_results{1}, contains all parameters used to run
% the function: "params" (this is a structure array itself).
% The second element in the cell array, bead_results{2}, contains the track
% results: "bead_final".
% "bead_final", is itself a structure array with end_frame by L elements,
% each of which is a structure with fields: CentreX, CentreY, ClipFlag, TooCloseToEdge, FrameNumber and SpotNumber.
% Along the dimension L, we have the different bead centres found on each
% frame (L will usually be the number of bead centres found in frame_start,
% or the largest number of spots ever accepted).
% For example, bead_final(100,8) is a structure with the above fields
% corresponding to the eighth found bead centre on frame 100.
%
% Example of params: bead_results{1}:
%
%             image_label: '25'
%             start_frame: 1
%               end_frame: 10
%      max_num_candidates: 200
%      subarray_halfwidth: 60
%     inner_circle_radius: 50
%                 rsq_min: 0.9000
%          d_coincid_cand: 3
%         d_coincid_found: 1
%                d_01_max: 30
%                d_02_max: 30
%                     rej: 10000
%               file_name: '210217r25.tif'
%               numFrames: 117
%             frame_Ysize: 1024
%             frame_Xsize: 1280
%
% bead_results1{2} is a struct array with fields:
%     CentreX
%     CentreY
%     ClipFlag
%     TooCloseToEdge
%     rsqFitX
%     rsqFitY
%     FrameNumber
%     BeadNumber
%     TrajNumber


%% DEFINITIONS and PARAMETERS:

% % One can save all parameter values in a file, e.g. "paramsForFindTrajects.m"
% % in the current directory. Just calling the name of that file loads the
% % parameter values into the workspace:
% paramsForFindTrajects
% % In this way, one can save different parameter sets for different data sets
% % in an easy manner and run two Matlabs at the same time working with different parameter sets.

% Alternatively, one can set PARAMETER values within this function below:  
%
% % Define data directory:
% % dir_data = 'Z:\Leake\Heiko Data\';
% % dir_data = 'C:\Isabel\ExperimData\HeikoData\';
% % Print data directory on command window to guide user:
% disp(' ') % empty line
% disp(['The data directory (.sif images) is: ',cd]) % display current directory.

% Save input parameters to "params" structure (for output):
params.image_label = image_label;
params.start_frame = start_frame;
params.end_frame = end_frame;

% Maximum number of candidate beads (if eg. 260000 candidate positions are
% found, we get an error in function pdist: "Distance matrix has more
% elements than the maximum allowed size in MATLAB"), hence, we limit the
% max number of candidates:
max_num_candidates = 200; % should be around 200.
% Save parameters to results as structure "params":
params.max_num_candidates = max_num_candidates;

% PARAMETERS for finding bead centres (see findBeadCentre1frame.m, inputs to the function):
subarray_halfwidth = 60; % (Default: 60 pixels). Halfwidth of image square subarray which includes bead and background around it.
inner_circle_radius = 50; % (Default: 50 pixels). Radius of circular mask around bead. 
% Save parameters to results as structure "params":
params.subarray_halfwidth = subarray_halfwidth;
params.inner_circle_radius = inner_circle_radius;

% PARAMETERS for deciding if we accept a bead centre found by the function
% findBeadCentre1frame or not:
rsq_min = 0.9; % minimum acceptable r-square value (0.9) 
% Goodness of parabolic fit to central peaks in cross-correlation of bead profile along x and y.
% Save parameters to results as structure "params":
params.rsq_min = rsq_min;

% % PARAMETER to decide if deflation method is applied or not (subtract each
% % found spot to allow detection of weaker spots):
% deflate = 1;
% % Save parameters to results as structure "params":
% params.deflate = deflate;

% PARAMETERS for eliminating coincident positions:
d_coincid_cand = 3; % distance (in pixels) for eliminating coincidences in bead-position candidates. Default: 3 pixels.
d_coincid_found = 1; % distance for eliminating coincidences in found bead centres.
% Save parameters to results as structure "params":
params.d_coincid_cand = d_coincid_cand; % distance for eliminating coincidences in spot candidates.
params.d_coincid_found = d_coincid_found; % distance for eliminating coincidences in found spot centres.

% PARAMETERS for building trajectories:
% For linking bead centres in current and previous frames:
d_01_max = 30; % max distance in pixels between bead centres in current and previous frames, for linking them into a trajectory (Default: 30 pix).
% Save parameters to results as structure "params":
params.d_01_max = d_01_max; 

% For linking loose bead centres in current frame and 2 frames ago (jump of 1 frame in trajectory):
d_02_max = 30; % max distance in pixels between bead centres in current frame and 2 frames ago. (default: 30 pix).
% Save parameters to results as structure "params":
params.d_02_max = d_02_max; 

% Use a very large number (larger than image size in pixels) for rejected asignments:
rej = 10000;
params.rej = rej; % Save parameters to results as structure "params".

% % Parameter to exclude a region from accepting spots (see later on, lines 322 and 507):
% exclude_region = 0;
% % Note that when exclude_region is 1, at the moment is it set for an image
% % with two channels, where spots are excluded also from around a horizontal
% % line centred on the image.
% exclude_region_width = subarray_halfwidth; % exclude edges of image, only look at spots with a centre far enough from image edges.
% params.exclude_region = exclude_region;
% params.exclude_region_width = exclude_region_width;

% Alternative way of selecting image sequence file:
% uigetfile opens a file dialog box to choose data file:
% [file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
% strcat('data (.sif image):','  ',path_data,file_data)
% open a file dialog box to choose analysis file:
% [file_analysis,path_analysis] = uigetfile({'*.xls'}, 'Chose analysis file (trajectory):');
% strcat('analysis file (.xls trajectory):','  ',path_analysis,file_analysis)
% 
% disp(' ') % empty line
% disp(['The start frame for finding bright spot trajectories will be ',num2str(start_frame)]) % start_frame is an input.
% disp(['The end frame for finding bright spot trajectories will be ',num2str(end_frame)]) % end_frame is an input.
% -----------------------------------------------------


%% Read in the image-sequence data:

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.
% --------------------------------------------------------------

% Save to parameters (for output):
params.file_name = image_path;
params.numFrames = numFrames;
params.frame_Ysize = frame_Ysize;
params.frame_Xsize = frame_Xsize;


%% Find candidate bead positions for start_frame (first frame), and find bead centres for those:

frame = image_data(start_frame).frame_data; % extract matrix data for first frame.
frame = double(frame);
disp(['frame number: ',num2str(start_frame)]) % print frame number to Command Window.

% Xpos is a matrix of the same size as frame, containing x values for all
% pixels and similarly for Ypos (used in future sections):
[Xpos,Ypos] = meshgrid(1:frame_Xsize,1:frame_Ysize);

[candidate_X_000,candidate_Y_000] = findCandidateBeadPositions(frame,1); % Second input: use method 1, which seems to work better.
% candidate_X_000 and candidate_Y_000 are two column vectors of the same
% length containing the x and y coordinates of the candidate bead positions found on the image.
% They contain integer numbers: coordinates or pixel numbers which give
% position on image plane.

% % REGION OF INTEREST: accept only beads within ROI (UNUSED).
% % A SignalMask can potentially be used to select a region of interest in the future. For now, it does nothing. 
% SignalMask = ones(size(frame,1),size(frame,2));
% % Reject candidate spots outside ROI (signal mask) (to speed up algorithm):
% candidate_X_00 = []; % initialise empty vectors before loop.
% candidate_Y_00 = []; 
% for nn = 1:length(candidate_X_000)
%    % Only use candidates for which there is a 1 in the SignalMask image:
%    if SignalMask(candidate_Y_000(nn),candidate_X_000(nn)) == 1
%        candidate_X_00 = [candidate_X_00; candidate_X_000(nn)];
%        candidate_Y_00 = [candidate_Y_00; candidate_Y_000(nn)];
%    end
% end
% disp(['no. of new candidate bead positions on start frame: ',num2str(length(candidate_X_00))])

% Without using a ROI:
candidate_X_00 = candidate_X_000;
candidate_Y_00 = candidate_Y_000;

disp(['no. of initial bead candidate positions: ',num2str(length(candidate_X_00))])

% Error control:
% Limit the max number of candidate spots (if eg. 260000 candidate spots are
% found, we will get an error in function pdist: "Distance matrix has more
% elements than the maximum allowed size in MATLAB").
% Select only the first max_num_candidates then.
if length(candidate_X_00) > max_num_candidates
    candidate_X_00 = candidate_X_00(1:max_num_candidates);
    candidate_Y_00 = candidate_Y_00(1:max_num_candidates);
    disp(['NOTE!! no. of candidate bead positions has been limited to ',num2str(max_num_candidates)])
end

% % Check graphically:
% imshow(frame,[]);
% hold on;
% plot(candidate_X_00,candidate_Y_00,'*');
% figure;

% Eliminate candidate positions closer to each other than d_coincid_cand (3 pixels):
[candidate_X_0,candidate_Y_0,pos_to_keep1] = eliminateCoincidentPositions(candidate_X_00,candidate_Y_00,d_coincid_cand);

disp(['no. of total candidate bead positions after eliminating coincidences: ',num2str(length(candidate_X_0))])

% Subtract background from entire frame (method 1 is faster) considering bead candidate positions:
frameNoBgnd = removeBgnd(frame,candidate_X_0,candidate_Y_0,inner_circle_radius,subarray_halfwidth,1);

% Find bead centres for first frame and decide if we accept them or not:
n =1; % Initialise index n (index for accepted bead centre positions):
frame_to_search = frameNoBgnd; % Initialise frame to search for bead centres.

% Find bead centres with sub-pixel precision:
for m = 1:size(candidate_X_0,1) % loop through all candidate bead positions.
    % Now find bead centre using function findBeadCentre1frame:
    % use candidate bead positions as initial estimates and then refine to find bead centre with sub-pixel precision.
   
    bead_result = findBeadCentre1frame(frame_to_search,candidate_X_0(m),candidate_Y_0(m),inner_circle_radius,subarray_halfwidth);
    bead_result.FrameNumber = start_frame; % Add new field containing frame number (time) to result structure.
    
    % Accepted bead centres:
    % if (bead_result.ClipFlag == 0 && ... % use this to exclude beads close to image edge.
    if (bead_result.ClipFlag < 2 && ... % use this to accept all beads;
            bead_result.TooCloseToEdge == 0 && ...
            bead_result.rsqFitX >= rsq_min && ...
            bead_result.rsqFitY >= rsq_min)
            % Only accept and save result of found bead centre if clipping flag = 0 and if values of rsquare of fits are acceptable.
            bead_result.BeadNumber = n; % Add new field containing bead number to result structure.
            bead_final(start_frame,n) = bead_result; % store "good" found bead centres.
            % This is also saved in the final result bead_final, structure array.
            % first index is for frame number, second index is for bead number.
                      
            n = n+1; % advance index n for accepted bead centres.
    end
end

% % display the number of accepted bead centres for this frame:
disp(['no. of accepted bead centres in first frame: ',num2str(n-1)])

% Convert results of found bead-centre positions to a useful form :
% Error control: if no bead centres were accepted:
if (n-1) == 0 
    found_bead_CentreX = [];
    found_bead_CentreY = [];
    % Need to create the whole bead_final structure with all its fields
    % here, just in case the number of accepted spots in the first frame is
    % zero, in order not to get error: "Subscripted assignment between
    % dissimilar structures".
    % Save empty bead (we need this, otherwise if in the last frame the no. of accepted beads is 0, there will be no result bead_final(end_frame,:) and the following functions will fail).
    bead_final(start_frame,n).CentreX = [];
    bead_final(start_frame,n).CentreY = [];
    bead_final(start_frame,n).ClipFlag = [];
    bead_final(start_frame,n).TooCloseToEdge = []; 
    bead_final(start_frame,n).rsqFitX = [];
    bead_final(start_frame,n).rsqFitY = [];
    bead_final(start_frame,n).FrameNumber = [];
    bead_final(start_frame,n).BeadNumber = []; 
    bead_final(start_frame,n).TrajNumber = [];
else
    found_bead_CentreX = [bead_final(start_frame,:).CentreX]'; % column vector with found CentreX positions of all candidate spots.
    found_bead_CentreY = [bead_final(start_frame,:).CentreY]'; % column vector with found CentreY positions of all candidate spots.
end

% Check graphically:
figure;
imshow(frame,[]);
hold on;
plot(found_bead_CentreX,found_bead_CentreY,'o','Color','g','MarkerSize',10) % plot accepted spot centres in green.
pause(0.1); % this pause is needed to give time for the plot to appear
hold off;
% -----------------------------------------------------------------------



%% Loop through selected frames:

tr =1; % initialise trajectory index.

for k = (start_frame+1):end_frame
    % to go through all frames do instead: for k = 1:length(sifData)
    
    frame = image_data(k).frame_data; % extract frame data which is stored in field 'frame_data'.
    frame = double(frame);
    
    imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
    hold on;
    
    disp(['frame number: ',num2str(k)]) % print frame number to Command Window.
    
    %-------------------------------------
    % Find new candidate bead positions for this frame:      
    [candidate_X_000,candidate_Y_000] = findCandidateBeadPositions(frame,1); % Second input: use method 1, which seems to work better.
    % the subindex "_000" in candidate_X_000 indicates newly found bead
    % candidates for the current frame. 
    
    % % Reject candidate spots outside ROI (SignalMask) (to speed up algorithm):
    % candidate_X_00 = []; % initialise empty vectors before loop.
    % candidate_Y_00 = [];
    % for nn = 1:length(candidate_X_000)
    %    % Only use candidates for which there is a 1 in the SignalMask image:
    %    if SignalMask(candidate_Y_000(nn),candidate_X_000(nn)) == 1
    %        candidate_X_00 = [candidate_X_00; candidate_X_000(nn)];
    %        candidate_Y_00 = [candidate_Y_00; candidate_Y_000(nn)];
    %    end
    % end
    % disp(['no. of new candidate bead positions on start frame: ',num2str(length(candidate_X_00))])
    %------------------------------------- 
    
    % Without using a ROI:
    candidate_X_00 = candidate_X_000;
    candidate_Y_00 = candidate_Y_000;
    
    % Error control:
    % Limit the max number of candidate bead positions (if eg. 260000 candidates are
    % found, we will get an error in function pdist: "Distance matrix has more
    % elements than the maximum allowed size in MATLAB").
    % Select only the first max_num_candidates then.
    if length(candidate_X_00) > max_num_candidates
        candidate_X_00 = candidate_X_00(1:max_num_candidates);
        candidate_Y_00 = candidate_Y_00(1:max_num_candidates);
        disp(['NOTE!! no. of candidate bead positions has been limited to ',num2str(max_num_candidates)])
    end
    
    % % Check graphically:
    % imshow(frame,[]);
    % hold on;
    % plot(candidate_X_0,candidate_Y_0,'*');
    % figure;
    
    disp(['no. of initial bead candidate positions: ',num2str(length(candidate_X_00))])
    
    % Eliminate candidate positions closer to each other than d_coincid_cand (3 pixels):
    [candidate_X_0,candidate_Y_0,pos_to_keep1] = eliminateCoincidentPositions(candidate_X_00,candidate_Y_00,d_coincid_cand);
    
    disp(['no. of total candidate bead positions after eliminating coincidences: ',num2str(length(candidate_X_0))])
    
    % Subtract background from entire frame (method 1 is faster) considering bead candidate positions:
    frameNoBgnd = removeBgnd(frame,candidate_X_0,candidate_Y_0,inner_circle_radius,subarray_halfwidth,1);

    % Find bead centres and decide if we accept them or not:
    n =1; % Initialise index n (index for accepted bead centre positions):
    frame_to_search = frameNoBgnd; % Initialise frame to search for bead centres.

    % Find bead centres with sub-pixel precision:
    for m = 1:size(candidate_X_0,1) % loop through all candidate bead positions.
        % Now find bead centre using function findBeadCentre1frame:
        % use candidate bead positions as initial estimates and then refine to find bead centre with sub-pixel precision.
        
        bead_result = findBeadCentre1frame(frame_to_search,candidate_X_0(m),candidate_Y_0(m),inner_circle_radius,subarray_halfwidth);              
        % index k is for frame number, index m is for spot number
        bead_result.FrameNumber = k; % Add new field containing frame number (time) to result structure.
        
        % Accepted bead centres:
        % if (bead_result.ClipFlag == 0 && ... % to exclude beads close to edge
        if (bead_result.ClipFlag < 2 && ... % use this to accept all beads;
                bead_result.TooCloseToEdge == 0 && ...
                bead_result.rsqFitX >= rsq_min && ...
                bead_result.rsqFitY >= rsq_min)
            % Only accept and save result of found bead centre if clipping flag = 0 and if values of rsquare of fits are acceptable.
            
            bead(k,n) = bead_result; % store accepted found bead centres in this preliminary result.
            % first index is for frame number, second index is for spot number.
            %    %--------------------------------------------------
            %    plot(bead(k,n).CentreX,bead(k,n).CentreY,'o','Color','r','MarkerSize',10); % Plot found bead centres in red
            %    %--------------------------------------------------
                      
            n = n+1; % advance index n for accepted bead centres.
        end
    end
    
    % % display the number of accepted bead centres for this frame:
    disp(['no. of accepted bead centres: ',num2str(n-1)])
    
    % % The following two lines are used together with the previous two
    % "plot" and "imshow" (commented off) lines:
    pause(0.1); % this pause is needed to give time for the plot to appear
    %    hold off;
    
    % Convert results of found bead-centre positions to a useful form:
    if (n-1) == 0 % error control: if no spots were accepted.
        found_bead_CentreX = [];
        found_bead_CentreY = [];
        bead_final(k,n).BeadNumber = []; % Save empty bead (we need this, otherwise if in the last frame the no. of accepted beads is 0, there will be no result bead_final(end_frame,:) and the following functions will fail).
    else
        found_bead_CentreX = [bead(k,:).CentreX]'; % column vector with found CentreX positions of all candidate bead centres.
        found_bead_CentreY = [bead(k,:).CentreY]'; % column vector with found CentreY positions of all candidate bead centres.
        
        %-------------------------------
        % This is not really necessary but I keep it:
        % Eliminate coincidences in result of last found bead centres for a given frame (for distance <1):
        [found_bead_CentreX found_bead_CentreY pos_final] = eliminateCoincidentPositions(found_bead_CentreX,found_bead_CentreY,d_coincid_found);
        % pos_final contains positions of selected, kept bead centres.
        %-------------------------------
        
        % Save final bead centres to variable bead_final:
        n=1; % index for final kept spot.
        for ii = 1:length(pos_final)
            mientras = bead(k,pos_final(ii)); % intermediate result.
            mientras.BeadNumber = n; % Add new field containing spot number to result structure.
            if k > (start_frame+1)
                mientras.TrajNumber = []; % to avoid error of disimilar structures
            end
            bead_final(k,n)=mientras; % final result structure of accepted spot centres.
            n = n+1;
        end
        % Plot found spot centres:
        pause(0.5);
        plot(found_bead_CentreX,found_bead_CentreY,'o','Color','g','MarkerSize',10) % plot final accepted spot centres in green.
        pause(0.1); % this pause is needed to give time for the plot to appear
        hold off;
        
        disp(['no. of final found bead centres after eliminating coincidences: ',num2str(length(found_bead_CentreX))])
    end
    
    
    
    %--------------------
    % LINKING SPOTS INTO TRAJECTORY SEGMENTS:
    
    % Link found and accepted spots into trajectory segments:
    
    % Trajectory index tr is initialised to 1 outside the loop through frames (k loop).
    
    % Do differently FOR SECOND FRAME (k == start_frame+1): compare only accepted spots in
    % previous and current frames:
    if k == start_frame+1 && ... % If second frame and
            (n-1)~=0 && ... % if the number of accepted spot centres is not zero and
            isempty([bead_final(k-1,:).BeadNumber])==0 && ... % at least 1 accepted spot in previous frame and
            isempty([bead_final(k,:).BeadNumber])==0 % at least 1 accepted spot in current frame.
        % There are no trajectories jet, so compare accepted spots in previous and current frames:
        N0 = max(cat(1,bead_final(k-1,:).BeadNumber));  % no. of accepted spots in previous frame.
        % Note: cat(1,bead_final(k-1,:).BeadNumber) gives a column vector with the values of SpotNumber for all non-empty accepted spots in frame k-1.
        N1 = max(cat(1,bead_final(k,:).BeadNumber));  % no. of accepted spots in current frame.
        
        % Create cell arrays with empty elements to pre-asign sizes:
        d01 = cell(N0,N1); % Note: d01 is a cell array (matrix) but d_01 below is a scalar.
        Iratio01 = cell(N0,N1); % Note: Iratio01 is a cell array (matrix) but Iratio_01 below is a scalar.
        SigmaRatio01 = cell(N0,N1); % Note: SigmaRatio01 is a cell array (matrix) but SigmaRatio_01 below is a scalar.
        
        for q0 = 1:N0 % loop though accepted bead centres in previous frame.
            for q1 = 1:N1 % loop though accepted bead centres in current frame.
                % d_01: distance between bead centres in previous and current frames:
                d_01 = sqrt((bead_final(k-1,q0).CentreX-bead_final(k,q1).CentreX)^2+(bead_final(k-1,q0).CentreY-bead_final(k,q1).CentreY)^2);
                %                 d_01
              
                % Accept and save trajectory if bead centres in previous and
                % current frames fulfill the following conditions:
                if d_01 < d_01_max   % see PARAMETERS at start of this function.           
                    % Asign accepted values to cell array elements to store them:
                    d01{q0,q1} = d_01; % use {} for cell arrays.
                else % rejected asignments:
                    d01{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                end
            end
            
%                         d01
%                         [d01{q0,:}]
            % Note that [d01{q0,:}] gives only non-empty elements of row q0
            % in the cell array d01 as a row vector, that's why we had to
            % give a numeric value rej to non-accepted asignments.
            
            % Note that if all asignments in previous step are rejected,
            % [d01{q0,:}] will be a list of rej values, and its minimum will
            % be rej.
            % If list of "linkable" beads, [d01{q0,:}], has no accepted
            % asignments (all values are rej):
            if min([d01{q0,:}]) == rej
                % Asign trajectory number 0 to the bead centre in the previous frame only:
                bead_final(k-1,q0).TrajNumber = 0;
            else % if there is at least one accepted asignment for a given bead centre in the previous frame:
                
                % Decide of all possible accepted bead centres (in current frame)
                % that could be linked to bead centre q0 in previous frame, which one is the best:
                % We take the best as the closest one to spot q0:
                q1_chosen = find([d01{q0,:}] == min([d01{q0,:}])); % find position of the minimum pair-wise distance.
                
%                             q1_chosen
                
                % Check if there is a better competing asignment for a given bead centre q1 in the current
                % frame from another bead centre in the previous frame.
                % Hence, check also column-wise in matrix d01 to avoid asigning a traj
                % number to a bead centre q1 in the current frame that had already
                % had a traj number asigned to it linking it to a different bead centre q0 in
                % the previous frame, which might be at a shorter distance
                % from it than the current one.
                
%                             [d01{:,q1_chosen}] % chosen column of d01 matrix of distances.
                
                % Asign trajectory numbers to structure bead_final:
                % If the found distance in that column is not the minimum one:
                if q0 ~= find([d01{:,q1_chosen}] == min([d01{:,q1_chosen}]));
                    bead_final(k-1,q0).TrajNumber = 0; % asign trajectory number 0 to bead centre in previous frame.
                else
                    bead_final(k-1,q0).TrajNumber = tr; % asign trajectory number to bead centre in previous frame, to bead_final structure.
                    bead_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to bead centre in current frame.
                    tr = tr+1; % advance trajectory-number index.
                end
            end
        end
        
        
        
    else % for FRAMES k >= start_frame+2, from third chosen frame on:
        
        % A) Compare loose bead centres (TrajNumber is 0) two frames ago (k-2)
        % to found bead centres in current frame (TrajNumber is []):
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % DO maybe!!
        if k ~= start_frame+1 && (n-1)~=0 && ... % If the number of accepted bead centres is not zero and
                isempty([bead_final(k-2,:).BeadNumber])==0 && ... % at least 1 accepted bead centre 2 frames ago and
                isempty([bead_final(k,:).BeadNumber])==0 % at least 1 accepted bead centre in current frame.
            
            N0 = max(cat(1,bead_final(k-2,:).BeadNumber));  % no. of accepted bead centres 2 frames ago.
            N1 = max(cat(1,bead_final(k,:).BeadNumber));  % no. of accepted bead centres in current frame.
            
            % Create cell arrays with empty elements to pre-asign sizes.
            d02 = cell(N0,N1); % Note: d02 is a cell array (matrix) but d_02 below is a scalar.
            
            for q0 = 1:N0 % loop though loose accepted bead centres 2 frames ago.
                if bead_final(k-2,q0).TrajNumber == 0 % only for loose (unlinked) bead centres (TrajNumber=0) two frames ago (so only rows q0 in matrix d01 which have unlinked bead centres will fill up).
                    for q1 = 1:N1 % loop though accepted bead centres in current frame.
                        % d_01: distance between bead centres in previous and current frames:
                        d_02 = sqrt((bead_final(k-2,q0).CentreX-bead_final(k,q1).CentreX)^2+(bead_final(k-2,q0).CentreY-bead_final(k,q1).CentreY)^2);
                        
                        %                         d_02
                        
                        % Accept and save trajectory if bead centres in previous and
                        % current frames fulfill the following conditions:
                        if d_02 < d_02_max   % see PARAMETERS at start of this function.                                
                            % Asign accepted values to cell array elements to store them:
                            d02{q0,q1} = d_02; % use {} for cell arrays.
                        else % rejected asignments:
                            d02{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                        end
                    end
                    
%                                         d02
%                                         [d02{q0,:}]
                    % Note that [d02{q0,:}] gives only non-empty elements of the cell array d01 as a row vector.
                    
                    % Note that if all asignments in previous step are rejected,
                    % [d02{q0,:}] will be a list of rej values, and its minimum will be rej.
                    % If list of "linkable" spots, [d02{q0,:}], has no accepted asignments (all values are rej):
                    if min([d02{q0,:}]) == rej
                        % Asign trajectory number 0 to the bead centre two frames ago only:
                        bead_final(k-2,q0).TrajNumber = 0; % point stays loose (unlinked).
                        
                    else % if there is at least one accepted asignment for a given bead centre two frames ago:
                        
                        % Decide of all possible accepted/saved bead centres (in current frame)
                        % that could be linked to bead centre q0 in previous frame (of all possible asignments), which one is the best:
                        % We take the best as the closest one to bead centre q0:
                        q1_chosen = find([d02{q0,:}] == min([d02{q0,:}])); % find position of the minimum pair-wise distance.
                        
%                                             q1_chosen
                        
                        % Check if there is a better competing asignment for a given bead centre q1 in the current
                        % frame from another bead centre q0 two frames ago.
                        % Hence, check also column-wise in d01 to avoid asigning a traj
                        % number to a bead centre q1 in the current frame that had already
                        % had a traj number asigned to it linking it to a different bead centre q0 two frames ago which might be at a shorter distance
                        % from it than the current one.
                        
%                                             [d02{:,q1_chosen}] % chosen column of d01 matrix of distances.
                        
                        % Asign trajectory numbers to structure bead_final:
                        % If the found distance in that column is not the minimum one:
                        if q0 ~= find([d02{:,q1_chosen}] == min([d02{:,q1_chosen}]));
                            % Asign trajectory number 0 to the bead centre two frames ago only:
                            bead_final(k-2,q0).TrajNumber = 0; % point stays loose (unlinked).
                        else
                            bead_final(k-2,q0).TrajNumber = tr; % asign trajectory number to bead centre two frames ago, to bead_final structure.
                            bead_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to bead centre in current frame.
                            tr = tr+1; % advance trajectory-number index.
                        end
                    end
                end
            end
        end
        
        
        % B) Compare loose bead centres (TrajNumber is []) and trajectories (TrajNumber is >0)
        % in previous frame (k-1) to found bead centres in current frame:
        if (n-1)~=0 && ... % If the number of accepted bead centres is not zero and
                isempty([bead_final(k-1,:).BeadNumber])==0 && ... % at least 1 accepted bead centre in previous frame.
                isempty([bead_final(k,:).BeadNumber])==0 % at least 1 accepted bead centre in current frame.
            
            % zzzzzzzzzzzzzzzzzzzzzzz
            N0 = max(cat(1,bead_final(k-1,:).BeadNumber));  % no. of accepted bead centres in previous frame.
            N1 = max(cat(1,bead_final(k,:).BeadNumber));  % no. of accepted bead centres in current frame.
            
            % Create cell arrays with empty elements to pre-asign sizes.
            d01 = cell(N0,N1);
            
            for q0 = 1:N0 % loop though accepted bead centres in previous frame.
                for q1 = 1:N1 % loop though accepted bead centres in current frame.
                    % d_01: distance between spot centres in previous and current frames:
                    d_01 = sqrt((bead_final(k-1,q0).CentreX-bead_final(k,q1).CentreX)^2+(bead_final(k-1,q0).CentreY-bead_final(k,q1).CentreY)^2);
                     
                    %                     d_01
                    
                    % Accept and save trajectory if bead centres in previous and
                    % current frames fulfill the following conditions:
                    if d_01 < d_01_max   % see PARAMETERS at start of this function.
                        % Asign accepted values to cell array elements to store them:
                        d01{q0,q1} = d_01; % use {} for cell arrays.
                    else % rejected asignments:
                        d01{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                    end
                end
                
%                                 d01
%                                 [d01{q0,:}] % last row of d01 matrix of distances.
                % Note that [d01{q0,:}] gives only non-empty elements
                % of row q0 in the cell array d01, as a row vector.
                
                % Note that if all asignments in previous step are rejected,
                % [d01{q0,:}] will be a list of rej values, and its minimum will
                % be rej.
                % If list of "linkable" spots, [d01{q0,:}], has no accepted
                % asignments (all values are rej):
                if min([d01{q0,:}]) == rej
                    
                    if isempty(bead_final(k-1,q0).TrajNumber) % if point in previous frame was not part of a trajectory (TrajNumber=[]):
                        % Asign trajectory number 0 to the bead centre in the previous frame only:
                        bead_final(k-1,q0).TrajNumber = 0;                       
                    end
                    
                else % if there is at least one accepted asignment for a given bead centre two frames ago:
                    
                    % Decide of all possible accepted/saved bead centres (in current frame)
                    % that could be linked to bead centre q0 in previous frame, which one is the best:
                    % We take the best as the closest one to spot q0:
                    q1_chosen = find([d01{q0,:}] == min([d01{q0,:}])); % find position of the minimum pair-wise distance.
                    
%                                     q1_chosen
                    
                    % Check if there is a better competing asignment for a given bead centre q1 in the current
                    % frame from another bead centre q0 in the previous frame.
                    % Hence, check also column-wise in d01 to avoid asigning a traj
                    % number to a bead centre q1 in the current frame that had already
                    % had a traj number asigned to it linking it to a different bead centre q0 in
                    % the previous frame which might be at a shorter distance
                    % from it than the current one.
                    
%                                     [d01{:,q1_chosen}]  % chosen column of d01 matrix of distances.
                    
                    % If the found distance in that column is not the minimum one:
                    if q0 ~= find([d01{:,q1_chosen}] == min([d01{:,q1_chosen}]));
                        if isempty(bead_final(k-1,q0).TrajNumber) % if point in previous frame was not part of a trajectory (TrajNumber=[]):
                            % Asign trajectory number 0 to the bead centre in the previous frame only:
                            bead_final(k-1,q0).TrajNumber = 0;
                        end
                    else
                        % Asign trajectory numbers to structure bead_final:
                        if bead_final(k-1,q0).TrajNumber > 0 % if point in previous frame was already part of a trajectory:
                            bead_final(k,q1_chosen).TrajNumber = bead_final(k-1,q0).TrajNumber; % asign that trajectory number to bead centre in current frame.
                        else % if point in previous frame was not part of a trajectory:
                            bead_final(k-1,q0).TrajNumber = tr; % asign new trajectory number to bead centre in previous frame.
                            bead_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to bead centre in current frame.
                            tr = tr+1; % advance trajectory-number index.
                        end
                    end
                end
            end
            % zzzzzzzzzzzzzzzzzzzzzzz
        end
    end
    %-----------------
    
    
end  % loop through selected frames



%% OUTPUT OF SPOT-FINDING PROCESS: final output bead_results:
%
bead_results = {params, bead_final};
% params is a structure array containing all parameters used to run the
% function.
% bead_final, is a structure array with end_frame x L elements,
% each of which is a structure with fields:
%     'CentreX'
%     'CentreY'
%     'rsqFit'
%     'ClipFlag'
%     'TooCloseToEdge' 
%     'TrajNumber'
%     'FrameNumber'
%     'SpotNumber'
%
% Along the dimension L, we have the different bead centres found on each
% frame (L will often be the number of bead centres found in frame_start,
% it is always the largest number of spots ever accepted on one frame).
% For example, bead_final(100,8) is a structure with the above fields
% corresponding to the eighth found bead centre on frame 100.
%
% Note that even if we only analyse from start_frame to end_frame,
% bead_final is a list containing empty structure arrays from index 1 to
% index start_frame, and then the found bead centres for the analysed frames start_frame to end_frame.
%
% The result is padded to a fixed number of bead centre structures for each
% frame (the maximum no. of accepted found bead centres of all frames), so that for a given
% frame in which less found bead centres have been accepted, the remaining
% elements are padded with empty structures with empty fields [].
% To check if a given bead centre is empty: isempty(bead_final(101,10).CentreX)
% gives 1 if field "CentreX" of tenth bead centre found and accepted in frame 101
% is empty (equal to []).
%
% e.g. cat(1,bead_final(100,:).SigmaFit) gives a vector column with all the
% non-empty SigmaFit values for all bead centres accepted in frame 100.
