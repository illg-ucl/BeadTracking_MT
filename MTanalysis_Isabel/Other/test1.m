function result = test1()
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
% Generate synthetic (x,y) data arranged as a circle in x-y plane and add
% noise; Then fit to a circle to find centre and radius of circle;
%
% INPUTS:
% 
arrayHalfSize = 50; % Half size of total array generated (to plot and detect circle);
xcentre = 51; % centre of circle along x;
ycentre = 51; % centre of circle along y;
radius = 20; % radius of circle;
Ndata = 50; % number of (x,y) points of data to generate along the circle;
noiseLevel = 1; % noise level amplitude factor.

% Pre-allocate size, column vector arrays:
circleXdata = zeros(Ndata,1);
circleYdata = zeros(Ndata,1);
% Go around circle to generate values in circle:
angleStep = (pi()/180)*360/Ndata; % angle step in rad.
% Loop to allocate values, including random noise (rand gives value between 0 and 1): 
for i = 1:Ndata
    circleXdata(i) = xcentre + radius*cos(angleStep*(i-1))+ noiseLevel*rand;
    circleYdata(i) = ycentre + radius*sin(angleStep*(i-1))+ noiseLevel*rand;
end

figure;
plot(circleXdata,circleYdata,'xr')
xlabel('x')
ylabel('y')
daspect([1 1 1]) % set aspect ratio to 1.

% Initialise matrix to plot circle:
mCircle = zeros(2*arrayHalfSize+1); 
% Fill in the matrix/image with ones at locations of circle points (rounded
% up):
for ii = 1:Ndata
    x_rounded = round(circleXdata(ii)); % first round to closest integer position;
    y_rounded = round(circleYdata(ii));
    mCircle(y_rounded,x_rounded) = 1; % assign value 1 to that coordinate.
end

figure;
imshow(mCircle,[]);
hold;
% Use Hough Transform to find circular edge and its centre and radius:
% note that it is rather sensitive to the 'Sensitivity' value, tweak until 
% it only finds one circle, it works well for 50 points in the circle using 0.97.
% Look in the range [0.7*radius, 1.3*radius]:
[centre_found, radius_found] = imfindcircles(mCircle,[round(0.7*radius) round(1.3*radius)],'Method','TwoStage','Sensitivity',0.97);
viscircles(centre_found, radius_found,'Color','r'); % display circle found with red circle overlaid;

% Error control:
if isempty(centre_found)
    disp('No circle found. Exiting function')
    %% Create partial output structure:
    result.xcentre = xcentre;
    result.ycentre = ycentre;
    result.radius = radius;
    result.circleXdata = circleXdata;
    result.circleYdata = circleYdata;
    result.xcentre_found = [];
    result.ycentre_found = [];
    result.radius_found = [];
else
    %% Create structure to store output results:
    result.xcentre = xcentre;
    result.ycentre = ycentre;
    result.radius = radius;
    result.circleXdata = circleXdata;
    result.circleYdata = circleYdata;
    result.xcentre_found = centre_found(2);
    result.ycentre_found = centre_found(1);
    result.radius_found = radius_found;

end




%% 
% % The following attempt to fit does not work:
% % Now, try to recover centre and radius from the two vectors
% % circleXdata,circleYdata by fitting to a circle.
% % The equation for the circle is (x-a)^2+(y-b)^2=r^2. So we can write it
% % also as y = b + sqrt(r^2-(x-a)^2) and fit y as a function of x, with fit
% % parameters r, a, b.
% fun_for_fit = fittype('sqrt(r1^2-(x-a1)^2)+b1','independent','x'); % define funtion to fit to, with 'x' as independent variable;
% options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% options.StartPoint = [radius+rand xcentre+rand ycentre+rand]; % give guess parameters for fit. This avoids a warning message. 
% options.Lower = [0.8*radius 0.8*xcentre 0.8*ycentre]; % Lower bounds for fit parameters (this is quite random).
% options.Upper = [1.2*radius 1.2*xcentre 1.2*ycentre]; % Upper bounds for fit parameters (this is quite random).
% options.MaxFunEvals = 3000;
% options.MaxIter = 3000;
% [fit_result gof] = fit(circleXdata,circleYdata,fun_for_fit,options); % do fit. 
% % fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the good of fitness.
% fit_param_names = coeffnames(fit_result); % fit parameter names to get their order.
% fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. 
% radius_fit = fit_param_values(1);
% xcentre_fit = fit_param_values(2); 
% ycentre_fit = fit_param_values(3); 
% rsq_fit = gof.rsquare; % rsquare coefficient of fit.
% errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
% errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).

