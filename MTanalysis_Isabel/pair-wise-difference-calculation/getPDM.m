function PDM = getPDM(x)
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
% Get pairwise difference matrix.
% Note that input vector x needs to be a column vector.

PDM = tril(repmat(x,[1 length(x)])-repmat(x',[length(x),1]),-1); 
% TRIL(X,K) is the elements on and below the K-th diagonal
% of X .  K = 0 is the main diagonal, K > 0 is above the
% main diagonal and K < 0 is below the main diagonal.
% REPMAT(A,[M N]) creates a large matrix consisting of an M-by-N
% tiling of copies of A. The size of B is [size(A,1)*M, size(A,2)*N].
end

% Examples for getPDM:
% x = [1 3 5 6]'; % column vector.
% repmat(x,[1 length(x)]) =
%      1     1     1     1
%      3     3     3     3
%      5     5     5     5
%      6     6     6     6
% repmat(x',[length(x),1]) =
%      1     3     5     6
%      1     3     5     6
%      1     3     5     6
%      1     3     5     6
% repmat(x,[1 length(x)])-repmat(x',[length(x),1]) =
%      0    -2    -4    -5
%      2     0    -2    -3
%      4     2     0    -1
%      5     3     1     0
% getPDM(x) = 
%      0     0     0     0
%      2     0     0     0
%      4     2     0     0
%      5     3     1     0
% The displacements over each time interval actually correspond to each
% diagonal in getPDM. For example, to extract first diagonal below the
% middle: 
% diag(getPDM(x),-1)= [2 2 1]';