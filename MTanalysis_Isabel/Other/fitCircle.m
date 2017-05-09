function result = fitCircle(xdata,ydata,radius_estimate,S)
% 
% % ========================================
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
% Fit (x,y) data arranged as a circle in x-y plane to find centre and radius;
%
% INPUTS:
% xdata: column vector with x positions;
% ydata: column vector with y positions;
% radius_estimate: estimated radius of circle in pixels.
% S: Sensitivity of the method for function imfindcircles. Use 0.97 and
% tweak around that value aimint to detect one circle in the image generated.

arrayHalfSize = 50; % Half size of total array generated (to plot and detect circle);
% Make sure arrayHalfSize is larger than the circle's radius.

plot(xdata,ydata,'xr')
xlabel('x')
ylabel('y')
daspect([1 1 1]) % set aspect ratio to 1.

% Initialise matrix to plot circle as image to then detect it:
mCircle = zeros(2*arrayHalfSize+1); 
% Fill in the matrix/image with ones at locations of circle points (rounded
% up):
for ii = 1:length(xdata)
    x_rounded = round(xdata(ii)); % first round to closest integer position;
    y_rounded = round(ydata(ii));
    mCircle(y_rounded,x_rounded) = 1; % assign value 1 to that coordinate.
end

figure;
imshow(mCircle,[]);
hold;
% Use Hough Transform to find circular edge and its centre and radius:
% note that it is rather sensitive to the 'Sensitivity' value, tweak until 
% it only finds one circle, it works well for 50-150 points in the circle using S=0.97.
% Look in the range [0.7*radius_estimate, 1.3*radius_estimate]:
[centre_found, radius_found] = imfindcircles(mCircle,[round(0.7*radius_estimate) round(1.3*radius_estimate)],'Method','TwoStage','Sensitivity',S)
viscircles(centre_found, radius_found,'Color','r'); % display circle found with red circle overlaid;

% Error control:
if isempty(centre_found)
    disp('No circle found. Exiting function')
    %% Create partial output structure:
    result.xcentre_found = [];
    result.ycentre_found = [];
    result.radius_found = [];
else
    %% Create structure to store output results:
    result.xcentre_found = centre_found(2);
    result.ycentre_found = centre_found(1);
    result.radius_found = radius_found;

end


