function scriptToRun
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

%% To do tests on a single frame:
frame1 = extract1frameB(1);

[x1,y1] = findCandidateBeadPositions(frame1,1);

% Eliminate candidate positions closer to each other than 3 pixels:
[new_x1,new_y1,pos_to_keep1] = eliminateCoincidentPositions(x1,y1,3);

% Subtract background, method 1 is faster:
frameNoBgnd1 = removeBgnd(frame1,new_x,new_y1,50,60,1);

% Test finding bead centre on single frame:
s1 = findBeadCentre1frame(frameNoBgnd1,new_x1(1),new_y1(1),50,60);
 

%% To analyse video sequences: 

% Find trajectories and output them to one excel file:
t25 = FindTrajectsBeads('25',1,10);
linkTrajSegmentsBeads('25',1,10,t25,'image25');

% Accept tracks manually after visual inspection:
good_tracks1 = goThroughBeadTracksVideo('25',1,'end',3); 