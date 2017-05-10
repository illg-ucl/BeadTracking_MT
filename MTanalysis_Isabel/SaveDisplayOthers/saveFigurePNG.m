function saveFigurePNG(folder_name_for_saving,figName)
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
% Save current figure as a .png file in a directory named "folder_name_for_saving",
% which already exists. The figure is saved with the entered name "figName".
% Both inputs are strings.
% It closes the figure window at the end.
%
% Example: figName = strcat('frameAvg',image_label), with image_label='554'.

cd(folder_name_for_saving); % change to directory in which results are saved.
% Export the current figure window at screen size as a png into folder
% folder_name_for_saving.
set(gcf, 'PaperPositionMode', 'auto')  % Use screen size. (gcf=get current figure)
% h = get(gcf);
print('-dpng','-r300',figName)  % add -r300 (to save at 300 dpi, higher resolution) after -dpng to control resolution.
cd('..'); % go back to previous directory.
close; % deletes the current figure (many open figures take too much memory and make Matlab crash).
    