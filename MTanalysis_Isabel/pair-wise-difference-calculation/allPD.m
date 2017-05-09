function PD_all = allPD(input_vector)
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
% Calculation of Pair-wise Differences (PD) for an input vector. 
% Note that all pair-wise differences are calculated for all pairs, not only for
% consecutive pairs of points.
%
% Pair wise differences are calculated from start to end of input_vector:
% positive/negative pair-wise differences corresponds to values
% increasing/decreasing, respectively, from start to end of input_vector.
%
% Inputs: 
% input_vector: needs to be a COLUMN vector.
%
% Outputs: 
% PD_all: column vector with all pair-wise differences.

PDM_all = getPDM(input_vector); % Pair-wise Differences Matrix. See getPDM.m.
% Extract all useful values, extract all diagonals from matrix:
PD_all = []; % initialise vector of pair-wise differences (PD).
for i=1:size(PDM_all,1)-1
    PD_all = [PD_all diag(PDM_all,-i)']; 
end
PD_all = PD_all'; % transpose into column vector.