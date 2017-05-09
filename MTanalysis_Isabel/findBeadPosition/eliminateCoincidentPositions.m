function [new_xx,new_yy,candid_pos_to_keep] = eliminateCoincidentPositions(xx,yy,limit_dist)
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
% In a list of (x,y) points (candidate bead positions), when two of them
% are closer than the input limit_dist, eliminate one of them.
% The function checks the distances between all pairs of points and removes
% those points (x,y) which are closer than limit_dist pixels
% (distance<limit_dist) to another point in the list.
%
% INPUTS: column vectors with the x and y coordinates of the points:
% - xx: columns of positions of candidate beads (horizontal). Column vector.
% - yy: rows of positions of candidate beads (vertical). Column vector.
% xx and yy must have the same length.
% - limit_dist: distance in pixels, see above.
% 
% OUTPUTS:
% - new_xx: column vector, list of x positions after eliminating coincidences.
% - new_yy: column vector, list of y positions after eliminating coincidences.
% - candid_pos_to_keep: list of positions kept within original input vectors.
% 


% Put together in a matrix with two columns the x and y positions of
% all candidates:
allCandidates = [xx yy];
P = size(allCandidates,1);

% Error control: need at least two points to compare:
if P<2
    new_xx = xx;
    new_yy = yy;
    if P==1
        candid_pos_to_keep = 1;
    elseif P==0
        candid_pos_to_keep = [];
    end
    return
end

distCandidates = pdist(allCandidates); % Calculate Euclidean distance between pairs of candidate points.
link_result = linkage(distCandidates); % Link candidate points according to their distances.
% link_result is a matrix where the first two columns identify the points
% that have been linked and the third column contains distance between them.

% Eliminate coincident candidates which are closer than limit_dist to each other: 
% Find pairs of candidates with distances above or equal to limit_dist and keep them all (see to_keep1). 
% Find pairs of candidates with distances below limit_dist and keep only one of the pair (see to_keep2).
% Find solitary spots with distances above or equal to limit_dist to any other spot (see to_keep3).
to_keep1 = find(link_result(:,3)>=limit_dist & link_result(:,1)<=P & link_result(:,2)<=P);
to_keep2 = find(link_result(:,3)<limit_dist & link_result(:,1)<=P & link_result(:,2)<=P);
to_keep3 = find(link_result(:,3)>=limit_dist & xor(link_result(:,1)<=P,link_result(:,2)<=P));
intermediate = link_result(to_keep3,[1 2]);

% For debugging:
% link_result(to_keep1,1)
% link_result(to_keep1,2)
% link_result(to_keep2,1)
% intermediate(intermediate<=P)

% candidates positions to be kept:
candid_pos_to_keep = sort([link_result(to_keep1,1); link_result(to_keep1,2); link_result(to_keep2,1); intermediate(intermediate<=P)]);

% Candidates we keep:
new_xx = xx(candid_pos_to_keep);
new_yy = yy(candid_pos_to_keep);

% disp(['no. of initial spots in list: ',num2str(length(xx))])
% disp(['no. of final spots after eliminating coincidences in list: ',num2str(length(new_xx))])