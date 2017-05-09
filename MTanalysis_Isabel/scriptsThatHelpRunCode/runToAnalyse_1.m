% runToAnalyse.m
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

%% Move into correct directory:

%  cd('C:\Isabel\DataAnalysis\Ollie's\cyoA_mCherry_OlliesData\October2012\power1p9mW')


%% Calculate frame averages:

% frame_avg = frameAverage(image_label,start_frame,end_frame,display_result,save_png)

% frameAverage('1ob',1,'end',1,1); % fussy, focus is not very good.
% frameAverage('2ob',1,'end',1,1); % fussy, focus is not very good.
% frameAverage('3ob',1,'end',1,1); % ok
% frameAverage('5ob',1,'end',1,1); % ~ ok
% frameAverage('6ob',1,'end',1,1); % ~ ok 
% frameAverage('7ob',1,'end',1,1); % ~ ok
% frameAverage('8ob',1,'end',1,1); % ok
% frameAverage('9ob',1,'end',1,1); % ok
% frameAverage('10ob',1,'end',1,1); % ok


%% Find trajectories and output them to one excel file:

% See paramsForFindTrajects.m.
% See paramsForLinkTrajSegments.m.

% s3 = FindTrajects('3ob',11,200);
% linkTrajSegments('3ob',11,200,s3,'cyoA_mCherry');
% 
% s5 = FindTrajects('5ob',11,200);
% linkTrajSegments('5ob',11,200,s5,'cyoA_mCherry');
% 
% s6 = FindTrajects('6ob',8,200);
% linkTrajSegments('6ob',8,200,s6,'cyoA_mCherry');
% 
% s7 = FindTrajects('7ob',12,200);
% linkTrajSegments('7ob',12,200,s7,'cyoA_mCherry');
% 
% s8 = FindTrajects('8ob',11,200);
% linkTrajSegments('8ob',11,200,s8,'cyoA_mCherry');
% 
% s9 = FindTrajects('9ob',10,200);
% linkTrajSegments('9ob',10,200,s9,'cyoA_mCherry');
% 
% s10 = FindTrajects('10ob',14,200);
% linkTrajSegments('10ob',14,200,s10,'cyoA_mCherry');
% 
% save 'resultStructures' 's*' % save all result structures in a .mat file.


%% Generate lists of good tracks after visual inspection:

% % goThroughTracksVideo(image_label,n_traj_start,n_traj_end,minPointsTraj)

% % This generates the file "good_track_nums_3ob.mat" inside the current
% directory, etc.:

% goThroughTracksVideo('3ob',1,'end',3) % 4 tracks only.
% goThroughTracksVideo('5ob',1,'end',3) % 1 track only.
% goThroughTracksVideo('6ob',1,'end',3) % 9 tracks only.
% goThroughTracksVideo('7ob',1,'end',3) % 0 long enough tracks.
% goThroughTracksVideo('8ob',1,'end',3) % 19 tracks only.
% goThroughTracksVideo('9ob',1,'end',3) %  0 long enough tracks.
% goThroughTracksVideo('10ob',1,'end',3) % 5 tracks only, only 2 tracks in the cell which does not wobble.

% Note: images '7ob' and '9ob' have no long enough tracks (with at least 3
% points).


%% Then use showManyTrajAnalysis2.m amd showTrajAnalysis2.m to produce one
% analysis excel file and graph per track:

% See paramsForShowTrajAnalysis2.m
% showManyTrajAnalysis2(image_label,n_traj_start,n_traj_end,start_frame,tsamp,pixelsize_nm,showVideo,minPointsTraj) 

% showManyTrajAnalysis2('3ob',1,'end',11,0.04,46,1,3)
% showManyTrajAnalysis2('5ob',1,'end',11,0.04,46,1,3)
% showManyTrajAnalysis2('6ob',1,'end',8,0.04,46,1,3)
% showManyTrajAnalysis2('8ob',1,'end',11,0.04,46,1,3)
% showManyTrajAnalysis2('10ob',1,'end',14,0.04,46,1,3)

%%
%  Now see "runToAnalyse_wholeCell_I.m"
% Also see
% C:\Isabel\DataAnalysis\Ollie's\cyoA_mCherry_OlliesData\24October2012\in vivo cyoa-mCherry\summaryAnalysis.doc
