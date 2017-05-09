function [xx_new yy_new] = movingAvg(xx,yy,window,plotYN)
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
% Moving average of a vector;
% 
% INPUTS: 
% xx and yy: column vectors;
% window: integer number, number of points to average over.
% plotYN: generate plots? Yes (1) or no (0).
%
% OUTPUTS:
% xx_new, yy_new: new averaged column vectors ;
% 

% put the two column vectors into a matrix to sort them:
Mtest = [xx yy];
MtestSorted = sortrows(Mtest); % sort by order of the first column;
xx_sorted = MtestSorted(:,1);
yy_sorted = MtestSorted(:,2);

xx_new = xx_sorted; 
yy_new = zeros(size(yy_sorted));
for i = 1:(length(yy_sorted)-(window-1))
    % moving average over a window of size "window":
    cummul = yy_sorted(i);
    for j = 1:(window-1)
        cummul = cummul + yy_sorted(i+j);    
    end
    yy_new(i) = cummul/window; % average
end

if plotYN == 1
    figure; plot(xx_sorted,yy_sorted,'.');
    hold;
    plot(xx_new,yy_new,'r.');
end
