function [new_xx,new_yy] = excludeRegions(xx,yy,list_xstart,list_xend,list_ystart,list_yend)
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
% Function to exclude certain regions of the image. Eliminates all xx, yy
% particle positions found that fall within certain regions. 
% The regions to exclude are given by start coordinates (x_start, y_start)
% and end coordinates (x_end, y_end) that delimit rectangular boxes on
% image.
%
% INPUTS: column vectors with the x and y coordinates of the points:
% - xx: columns of positions of candidate particles (horizontal). Column vector.
% - yy: rows of positions of candidate particles (vertical). Column vector.
% xx and yy must have the same length.
% - list_xstart: vector of starting horizontal x positions of regions to
% exclude;
% - list_xend: vector of end horizontal x positions of regions to
% exclude;
% - list_ystart: vector of starting vertical y positions of regions to
% exclude;
% - list_yend: vector of end vertical y positions of regions to
% exclude;
% Note that all list_xstart,list_xend,list_ystart and list_yend must be vectors of
% the same length.
%
% 
% OUTPUTS:
% - new_xx: column vector, list of x positions after eliminating positions within excluded regions.
% - new_yy: column vector, list of y positions after eliminating positions within excluded regions.
% 

% Get number of excluded regions: 
Nregions = length(list_xstart);

% Initialise positions of candidates to be kept:
candid_pos_to_keep = [];

% Loop through all input positions:
for i = 1:length(xx)
    for k = 1:Nregions
        if xx(i)>= list_xstart(k) && xx(i)<= list_xend(k) && yy(i)>= list_ystart(k) && yy(i)<= list_yend(k)
            eliminate(k) = 1;
        else
            eliminate(k) = 0;
        end
    end
    % if input position is not in exclusion regions, keep it:
    if sum(eliminate) == 0
        candid_pos_to_keep = [candid_pos_to_keep i]; % append position in input vector to be kept
    end
end

% Output, candidates to keep:
new_xx = xx(candid_pos_to_keep);
new_yy = yy(candid_pos_to_keep);