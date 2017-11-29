function bead = findBeadCentre1frame(frameNoBgnd,x_estimate,y_estimate,inner_radius,subarray_halfsize)
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
% Function to iteratively find the centre of a bead on an image with
% sub-pixel precision.
%
% INPUTS:
% - frameNoBgnd is a matrix containing the background-subtracted image data.
% Good background subtraction is crucial for the cross-correlation later
% on, since vectors are padded with zeros and jumps in intensity lead
% to incorrect cross-correlation.
% frameNoBgnd can be obtained doing, e.g.:
% frame1 = extract1frameB(1);
% [x1,y1] = findCandidateBeadPositions(frame1,1);
% frameNoBgnd1 = removeBgnd(frame1,x1,y1,50,60,2);

% - x_estimate, y_estimate are the initial estimated centre positions fed to the
% iterative method. These can be the output of function
% findBeadCentre1frame.
%
% - inner_radius: radius of inner circular mask in pixels, used to carry
% out a second background subtraction if necessary (usually not). Make sure this is smaller than
% subarray_halfsize.
%
% - subarray_halfsize: the algorithm selects this number of pixels above and below and to
% left and right of estimated centroid pixel to form the selected image
% subarray (I). This depends on the bead size, must be larger than bead size, default~35-70.
%
%
% OUTPUTS: values resulting from the iterative method. The output is a
% structure "bead" with the following fields (see end of file):
% - CentreX % x-centre result found.
% - CentreY  % y-centre result found.
% - rsqFitX % r squared of parabolic fit of central peak of x cross-correlation.
% - rsqFitY % r squared of parabolic fit of central peak of y cross-correlation.
% - ClipFlag % 1 if candidate was closer to edge of image than inner_radius.
% - TooCloseToEdge % 1 if candidate position closer to edge than d_min (see PARAMETERS below).
%
% Example of how to call this function: s1 = findBeadCentre1frame(frameNoBgnd1,271,364,50,60),
% where frame1 is the image array (single-frame data). Then to obtain results
% from the structure s1, do, for example: s1.CentreX, s1.CentreY, etc.

%-----------------
% PARAMETERS:
d_min = 10; % minimum subarray halfsize in pixels. If candidate bead centre too close to edge, return empty result.
% Choose method 1 or method 2 for averaging radial profiles.
% Method 1: averages radial profiles over hemispheres.
% Method 2: averages radial profiles over quarters, -45deg to 45deg, etc.
averagingMethod = 2; % Default = 2.
% Step size for rho and angle, for converting image to a polar grid:
angle_step = 2; % angle step in degrees. Default = 2.
radial_step = 1; % rho step in pixels. Default = 1.
% Additional background subtraction: 1 for yes, 0 for no.
additional_bgnd_subtract = 0; % Default = 0 if background-subtracted image is used as input.
% N_tofit: Half-number of points to fit to a parabola in peak of cross-correlations.
% Default = 3.
N_tofit = 3; % number of points to fit at either side of the centre in the central peak of the cross-correlation.
guessA_fitXY = 1; % guess for parameter A for parabolic fit of central peak of cross-correlation, for both X and Y fits
% Default guessA_fitXY for magnetic beads (Crick): 1.
guessB_fitXY = -0.03; % guess for parameter B for parabolic fit of central peak of cross-correlation, for both X and Y fits
% Default guessB_fitXY for magnetic beads (Crick): -0.005.
%-----------------

% % Aid plot: plot original frame and estimate centre position
% % (x_estimate,y_estimate) overlayed as red circle on top:
% figure;
% imshow(frameNoBgnd,[],'InitialMagnification',150);
% hold;
% plot(x_estimate,y_estimate,'o','Color','r','MarkerSize',12);
% title('current bead with centre estimate')

d = subarray_halfsize; % subarray halfsize. d needs to be at least equal to the radius of the inner mask, inner_radius. Default: d=8, inner_radius=5.


%% Error control:
if subarray_halfsize < inner_radius
    error('findBeadCentre1frame:one','check input parameters: "subarray_halfsize" cannot be smaller than "inner_radius"')
end

% If the bead candidate is at the edge of the image, flag it with a clipping flag and take a subarray around
% it with a smaller size, as large as possible until the edge is reached, so reasign d:
clipping_flag = 0;
tooCloseToEdge = 0;
d_top = subarray_halfsize;
d_bottom = subarray_halfsize;
d_left = subarray_halfsize;
d_right = subarray_halfsize;
% If the bead is at edge of image, take a smaller subarray around
% it but as large as possible until the edge is reached,
% so reasign the d values:
if (round(y_estimate)-d_top)<1
    d_top = round(y_estimate) - 1;
end
if (round(y_estimate)+d_bottom)>size(frameNoBgnd,1)
    d_bottom = size(frameNoBgnd,1)-round(y_estimate);
end
if (round(x_estimate)-d_left)<1
    d_left = round(x_estimate) - 1;
end
if (round(x_estimate)+d_right)>size(frameNoBgnd,2)
    d_right = size(frameNoBgnd,2)-round(x_estimate);
end
% Chose the minimum distance:
d = min([d_top d_bottom d_left d_right]);

if d < inner_radius
    clipping_flag = 1;
end

if d < d_min
    % Error control: return empty result structure as outuput if bead-centre candidate is closer
    % to edge than d_min:
    tooCloseToEdge = 1;
    
    bead.CentreX = []; % x-centre result found.
    bead.CentreY = []; % y-centre result found.
    bead.rsqFitX = []; % r squared of parabolic fit of central peak of x cross-correlation.
    bead.rsqFitY = []; % r squared of parabolic fit of central peak of y cross-correlation.
    bead.ClipFlag = clipping_flag; % 1 if candidate was closer to edge of image than inner_radius.
    bead.TooCloseToEdge = tooCloseToEdge; % 1 if bead candidate was closer to edge of image than d_min.
    
else
         
    %% Create image subarray (I) around bead:
    % Create image subarray of size (2*d+1)x(2*d+1) centered on (x_estimate,y_estimate) centroid estimate:
    I = frameNoBgnd(round(y_estimate)-d:round(y_estimate)+d,round(x_estimate)-d:round(x_estimate)+d);
    % % Aid plot:
    % figure;
    % imshow(I,[])
    
    % Create matrices containing the x and y positions corresponding to the
    % subarray I. Ys is a matrix of the same size as I, containing y values for
    % all pixels and similarly for Xs. In Ys, y varies along first dimension
    % (vertical) in matrix. In Xx, x varies along second dimension (horiz).
    [Ys,Xs] = ndgrid(round(y_estimate)-d:round(y_estimate)+d,round(x_estimate)-d:round(x_estimate)+d);
    % ndgrid is used because the interpolation function interpn later on only accepts ndgrid format.
    
    
    %% Additional background subtraction:
    if additional_bgnd_subtract == 1
        % Generate inner circular mask for subarray, to calculate average
        % background per pixel and subtract it:
        inner_mask = zeros(2*d+1);
        % Assign value 1 to pixels for which radius <= inner_radius. C = hypot(A,B) returns SQRT(ABS(A).^2+ABS(B).^2).
        inner_mask = inner_mask | hypot(Xs-x_estimate, Ys-y_estimate) <= inner_radius;
        % Background mask is the negative of the inner circle mask, i.e.,
        % not-inner_mask. It has zeros in centre and ones in surrounding pixels:
        bgnd_mask = double(~inner_mask); % make of class double (not binary) for what follows.
        Ibgnd = I.*bgnd_mask; % bgnd-only image.
        pos_bgnd = find(bgnd_mask==1); % positions of bgnd-only intensities in bgnd mask.
        % Get mean (or median) background intensity per pixel in bgnd region, to
        % exclude hot pixels or nearby bright beads, Ibg_avg:
        % Ibg_avg = median(Ibgnd(pos_bgnd)); % use median.
        Ibg_avg = mean(Ibgnd(pos_bgnd)); % use mean.
        % bg_noise_std = std(Ibgnd(pos_bgnd)); % standard deviation of matrix elements in bgnd region.
        % Calculate background-corrected subarray image:
        I2 = I-Ibg_avg;
    else
        % If additional background subtraction not necessary (typically) do:
        I2 = I;
    end
    
    %% Generate image in polar coordinates:
    % Generate polar coordinate grid of (r,theta) values. Then calculate the
    % corresponding (x,y) cartesian coordinates for those values and
    % interpolate the image subarray I to find the corresponding values at
    % coordinates (r,theta). This then allows the calculation of radial average
    % profiles for each bead quarter (top, bottom, left, right).
    % Procedure is similar (though not equal) to paper: "Non-bias-limited tracking of spherical
    % particles, enabling nanometer resolution at low magnification", M. T. J.
    % van Loenhout, ..., Cees Dekker, et al. Biophys. J. 102, 2362 (2012).
    %
    % Polar grid.
    % r values between 0 and d, in steps of 1 pixel:
    % theta values between 0 and 360 degrees, in steps of 1 degree [for
    % (x,y)=(d,1), the minimum angle is ~1 degree.
    % angle_step is the angle step in degrees. See PARAMETERS section above.
    % radial_step is the rho step in pixels. See PARAMETERS section above.
    [rho_grid,theta_grid] = ndgrid(0:radial_step:d,(0:angle_step:360)*pi()/180);
    % theta changes along horizontal dimension (2) of matrices, rho along vertical one (1).
    
    % Calculate the corresponding cartesian coordinates:
    % Note that x_grid, y_grid are matrices with non-integer values centered around (0,0):
    [x_grid,y_grid] = pol2cart(theta_grid,rho_grid); % pol2cart transforms polar coordinates into cartesian ones.
    % Now since cartesian coordinate axes are defined as: y vertical (pointing
    % down) and x horizontal (pointing right), note that positive y values
    % correspond to bottom half of bead (for theta 0 to 180deg) and negative y
    % values correspond to top half of bead (for theta 180 to 360deg).
    
    % Now interpolate (2D cubic interpolation) to find image intensity value at
    % the x_grid, y_grid positions that correspond to the polar grid.
    % Having the integer-valued grids Xs, Ys (their values are
    % positions in the overall array frame) we can subtract the bead-centre
    % estimate to centre the positions at zero, so Xs-round(x_estimate) and
    % Ys-round(y_estimate). For these centered X, Y grids, the corresponding
    % image intensity is I2. For the new, non-integer valued arrays x_grid,y_grid
    % (corresponding to the polar grid), we find out the new image intensities,
    % Ipolar via interpolation.
    % Image subarray in polar coordinates, centred at (0,0):
    Ipolar = interpn(Ys-round(y_estimate),Xs-round(x_estimate),I2,y_grid,x_grid,'cubic');
    % Dimensions of Ipolar are same as dimensions of rho_grid and theta_grid.
    
    % % Aid plots:
    % % Plots of the polar grid used:
    % figure;
    % polarplot(theta_grid,rho_grid,'Marker','o','MarkerSize',2);
    % % Plot overlaying the original cartesian grid and the new interpolated
    % % points that derive from the polar grid:
    % figure;
    % mesh(Xs-round(x_estimate),Ys-round(y_estimate),I2); % original cartesian grid.
    % hold on;
    % mesh(x_grid,y_grid,Ipolar,'LineStyle','None','Marker','o','MarkerSize',3,'MarkerFaceColor','Flat')
    % title('Interpolated polar grid points on original cartesian grid')
    % hold off;
    % % surface plot of new interpolated points that derive from the polar grid:
    % figure; surf(x_grid,y_grid,Ipolar);
    % title('Interpolated polar grid')
    %
    % % original image of bead (after background subtraction) as contour plot:
    % figure;
    % min0 = min(min(I2));
    % max0 = max(max(I2));
    % step0 = round(max0-min0)/200; % 200 shades of grey (can use 1000 but slow).
    % contourf(Xs-round(x_estimate),Ys-round(y_estimate),I2,'LineStyle','None','LevelList',(min0:step0:max0));
    % colormap gray;
    % title('Original image of bead')
    %
    % % Interpolated (re-sampled) polar image:
    % figure;
    % contourf(x_grid,y_grid,Ipolar,'LineStyle','None','LevelList',(min0:step0:max0));
    % colormap gray
    % title('Interpolated (re-sampled) polar image of bead')
    
    
    %% Average radial intensity profiles along vertical and horizontal directions:
    
    % % Overall average radial intensity profile (averaged across all angles in the polar
    % % image):
    % profileRad = sum(Ipolar,2)/size(Ipolar,2);
    % % Plot of average radial profile:
    % figure;
    % plot(profileRad)
    % xlabel('r (pixels)')
    % xlabel('Intensity (arb)')
    % title('Radial intensity profile of bead (average over all angles)')
    
    % Choose method 1 or method 2 for averaging radial profiles. See
    % averagingMethod in PARAMETERS section above.
    % Method 1: averages radial profiles over hemispheres.
    % Method 2: averages radial profiles over quarters, -45deg to 45deg, etc.
    
    if averagingMethod == 1
        % METHOD 1: exactly as in paper "Non-bias-limited tracking of spherical
        % particles, enabling nanometer resolution at low magnification", M. T. J.
        % van Loenhout, ..., Cees Dekker, et al. Biophys. J. 102, 2362 (2012).
        %
        % Average radial intensity profiles over each 90-degree quarter:
        % top-left(TL), top-right (TR), bottom-left (BL), bottom-right (BR).
        % Note that y axis points down, hence theta angle grows clockwise. The
        % angles that separate the quarters are:
        pos_angle1 = 1+round(0/angle_step); % position of angle in polar image matrix.
        pos_angle2 = 1+round(90/angle_step);
        pos_angle3 = 1+round(180/angle_step);
        pos_angle4 = 1+round(270/angle_step);
        pos_angle5 = 1+round(360/angle_step);
        % Separate quarters:
        I_BR = Ipolar(:,pos_angle1:pos_angle2); % bottom right
        I_BL = Ipolar(:,pos_angle2:pos_angle3); % bottom left
        I_TL = Ipolar(:,pos_angle3:pos_angle4); % top left
        I_TR = Ipolar(:,pos_angle4:pos_angle5); % top right.
        % Calculate average radial profile for each quarter:
        profile_BR = sum(I_BR,2)/size(I_BR,2);
        profile_BL = sum(I_BL,2)/size(I_BL,2);
        profile_TL = sum(I_TL,2)/size(I_TL,2);
        profile_TR = sum(I_TR,2)/size(I_TR,2);
        % Obtain right profile as average of TR and BR profiles:
        profile_R = (profile_BR + profile_TR)/2;
        % Obtain left profile as average of TL and BL profiles:
        profile_L = (profile_BL + profile_TL)/2;
        % And similarly for top and bottom profiles:
        profile_T = (profile_TL + profile_TR)/2;
        profile_B = (profile_BL + profile_BR)/2;
        % HORIZONTAL (X) profile: concatenate right and left profiles:
        profile_H = [flip(profile_L(2:length(profile_L))); profile_R];
        % VERTICAL (Y) profile: concatenate top and bottom profiles:
        profile_V = [flip(profile_T(2:length(profile_T))); profile_B];
        % radial distance vector common for both profile_H and profile_V:
        radial_distance = [-flip(rho_grid(2:size(rho_grid,1),1)); rho_grid(:,1)];
    end
    
    if averagingMethod == 2
        % METHOD 2 for producing vertical/horizontal radial profile averages:
        % Average VERTICAL radial intensity profile: average over 90-degree vertical quarters:
        % average over 45deg to 135deg (bottom quarter) and concatenate to average over 225deg
        % to 315deg (top quarter) in the polar image.
        % Note that y axis points down, hence theta angle grows clock wise.
        pos_angleV1 = 1+round(45/angle_step); % position of angle in polar image matrix.
        pos_angleV2 = 1+round(135/angle_step);
        pos_angleV3 = 1+round(225/angle_step);
        pos_angleV4 = 1+round(315/angle_step);
        % Separate top and bottom quarters of image:
        Ibottom = Ipolar(:,pos_angleV1:pos_angleV2);
        Itop = Ipolar(:,pos_angleV3:pos_angleV4);
        % Get averaged radial top and bottom intensity profiles:
        profile_bottom = sum(Ibottom,2)/size(Ibottom,2);
        profile_top = sum(Itop,2)/size(Itop,2);
        % Concatenate top and bottom averaged radial profiles to get the average
        % profile along the VERTICAL direction:
        profile_V = [flip(profile_top(2:length(profile_top))); profile_bottom];
        radial_distance = [-flip(rho_grid(2:size(rho_grid,1),1)); rho_grid(:,1)];
        
        % Average HORIZONTAL radial intensity profile: average over 90-degree horizontal quarters:
        % Average over 0 to 45deg together with 315 to 360deg (RHS quarter),
        % and concatenate to average over 135 deg to 225deg (LHS quarter) in the polar image:
        pos_angleH1_a = 1+round(0/angle_step); % position of angle in polar image matrix.
        pos_angleH2_a = 1+round(45/angle_step);
        pos_angleH1_b = 1+round(315/angle_step); % position of angle in polar image matrix.
        pos_angleH2_b = 1+round(360/angle_step);
        pos_angleH3 = 1+round(135/angle_step);
        pos_angleH4 = 1+round(225/angle_step);
        % Separate right and left quarters of image:
        Iright_a = Ipolar(:,pos_angleH1_a:pos_angleH2_a);
        Iright_b = Ipolar(:,pos_angleH1_b:pos_angleH2_b);
        Iright = [Iright_a Iright_b];
        Ileft = Ipolar(:,pos_angleH3:pos_angleH4);
        % Get averaged radial right and left intensity profiles:
        profile_right = sum(Iright,2)/size(Iright,2);
        profile_left = sum(Ileft,2)/size(Ileft,2);
        % Concatenate right and left averaged radial profiles to get the average
        % profile along the HORIZONTAL direction:
        profile_H = [flip(profile_left(2:length(profile_left))); profile_right];
    end
    
    
    %% Cross-correlation to find bead centre with sub-pixel precision:
    % Calculate cross correlation between bead intensity profile and its mirror image.
    % This peaks at twice the displacement of the bead from the subarray centre.
    % Do for vertical and horizontal averaged bead intensity profiles.
    
    % Vertical (Y) cross-correlation of profile_V and its mirror image,
    % flip(profile_V):
    max_lag = length(profile_H); % Maximum lag, specified as an integer scalar.
    % How far away to go on displacement in each direction to calculate correlation, in pixels.
    [crossCorr_Y,steps_corY] = xcorr(profile_V,flip(profile_V),max_lag,'coeff');
    % crossCorr_Y is the cross-correlation vector of the two input signals.
    % steps_corY is a vector with the lags at which the correlations are computed.
    % The normalisation option 'coeff' normalizes so that the autocorrelations
    % at zero lag equal 1: Rfg,coeff(m) = Rfg(m)/sqrt(Rff(0)*Rgg(0)), ie.,
    % divides by autocorrelation.
    
    % Horizontal (X) cross-correlation of profile_H and its mirror image,
    % flip(profile_H):
    % max_lag is same as before.
    [crossCorr_X,steps_corX] = xcorr(profile_H,flip(profile_H),max_lag,'coeff');
    
    % % Aid plots:
    % % plot of vertical average profile:
%     figure;
%     subplot(2,2,1)
%     plot(radial_distance,profile_V)
%     xlabel('y (pixels)')
%     ylabel('Intensity (arb)')
%     title('Vertical intensity profile of bead (average over 90deg vertical cones)')
%     % plot of vertical average profile:
%     subplot(2,2,2)
%     plot(radial_distance,profile_H)
%     xlabel('x (pixels)')
%     ylabel('Intensity (arb)')
%     title('Horizontal intensity profile of bead (average over 90deg horizontal cones)')
%     % plot of vertical cross-correlation:
%     subplot(2,2,3)
%     plot(steps_corY,crossCorr_Y,'x')
%     xlabel('y (pixels)')
%     ylabel('normalised cross correlation (arb)')
%     % plot of horizontal cross-correlation:
%     subplot(2,2,4)
%     plot(steps_corX,crossCorr_X,'x')
%     xlabel('x (pixels)')
%     ylabel('normalised cross correlation (arb)')
    
    
    %% Use obtained cross-correlations to get correction to bead-centre position:
    % Along horizontal and vertical directions, with sub-pixel precision.
    
    % Fit cross-correlation to a parabola along X:
    % N_tofit is the number of points to fit at either side of the centre in
    % the central peak of the cross-correlation. Default is 3. See PARAMETERS
    % section at start of function.
    mid_pos = 1+(length(crossCorr_X)-1)/2; % position of middle point in cross-correlation.
    x_data = steps_corX(mid_pos-N_tofit:mid_pos+N_tofit)'; % data to fit, independent variable.
    y_data = crossCorr_X(mid_pos-N_tofit:mid_pos+N_tofit); % data to fit.
    fun_for_fit = fittype('A+B*(x-x0)^2','independent','x'); % define parabolic funtion to fit to, with 'x' as independent variable;
    options = fitoptions('Method','NonlinearLeastSquares'); % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
    options.StartPoint = [guessA_fitXY guessB_fitXY 0]; % give guess parameters for fit. This avoids a warning message.
    % guess for I0 is the max value of I in the background corrected image subarray.
    % options.Lower = [min_A min_B min_x0]; % Lower bounds for fit parameters.
    % options.Upper = [max_A max_B max_x0]; % Upper bounds for fit parameters.
    options.MaxFunEvals = 3000;
    options.MaxIter = 3000;
    [fit_result gof] = fit(x_data,y_data,fun_for_fit,options); % do fit. fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the good of fitness.
    % fit_param_names = coeffnames(fit_result); % fit parameter names, to get
    % their order: A, B, x0.
    fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit.
    A_fitX = fit_param_values(1);
    B_fitX = fit_param_values(2);
    x0_fitX = fit_param_values(3);
    rsq_fitX = gof.rsquare; % rsquare coefficient of fit.
    % errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    % errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    
    % % Aid plots:
%     x_data2 = min(x_data):(max(x_data)-min(x_data))/50:max(x_data); % plot fitted curve with higher sampling.
%     y_fitted = A_fitX+B_fitX*(x_data2-x0_fitX).^2;
%     figure;
%     subplot(1,2,1)
%     plot(x_data,y_data,'o');
%     hold on;
%     xlabel('x corr')
%     ylabel('cross-correl along X')
%     plot(x_data2,y_fitted);
%     title(strcat('x0 = ',num2str(x0_fitX),' pix'))
    
    % Fit cross-correlation to a parabola along Y:
    % N_tofit same as before.
    mid_pos = 1+(length(crossCorr_Y)-1)/2; % position of middle point in cross-correlation.
    x_data = steps_corY(mid_pos-N_tofit:mid_pos+N_tofit)'; % data to fit, independent variable.
    y_data = crossCorr_Y(mid_pos-N_tofit:mid_pos+N_tofit); % data to fit.
    % fun_for_fit and options same as before: fun_for_fit = fittype('A+B*(x-x0)^2','independent','x');
    [fit_result gof] = fit(x_data,y_data,fun_for_fit,options); % do fit. fit_result contains the fit coefficient values and their confidence intervals and "gof" gives the good of fitness.
    fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit.
    A_fitY = fit_param_values(1);
    B_fitY = fit_param_values(2);
    x0_fitY = fit_param_values(3);
    rsq_fitY = gof.rsquare; % rsquare coefficient of fit.
    % errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    % errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    
    % % Auxiliary plots:
    %     x_data2 = min(x_data):(max(x_data)-min(x_data))/50:max(x_data); % plot fitted curve with higher sampling.
    %     y_fitted = A_fitY+B_fitY*(x_data2-x0_fitY).^2;
    %     subplot(1,2,2)
    %     plot(x_data,y_data,'o');
    %     hold on;
    %     xlabel('y corr')
    %     ylabel('cross-correl along Y')
    %     plot(x_data2,y_fitted);
    %     title(strcat('y0 = ',num2str(x0_fitY),' pix'))

    % The displacements x0_fitX and x0_fitY obtained from fitting the
    % cross-correlation peaks are equal to twice the actual displacement of the
    % bead centre from the centre of the subarray:
    delta_X = x0_fitX/2;
    delta_Y = x0_fitY/2;
    
    x_new = round(x_estimate)+delta_X;
    y_new = round(y_estimate)+delta_Y;
    
    % % Return result:
%     disp(['new x_centre estimate = ',num2str(x_new)])
%     disp(['new y_centre estimate = ',num2str(y_new)])
    
    
    %% Output of the function:
    % The output is a structure "bead" with the following fields:
    bead.CentreX = x_new; % x-centre result found.
    bead.CentreY = y_new; % y-centre result found.
    bead.rsqFitX = rsq_fitX; % r squared of parabolic fit of central peak of x cross-correlation.
    bead.rsqFitY = rsq_fitY; % r squared of parabolic fit of central peak of y cross-correlation.
    bead.ClipFlag = clipping_flag; % 1 if candidate was closer to edge of image than inner_radius.
    bead.TooCloseToEdge = tooCloseToEdge; % 1 if bead candidate was closer to edge of image than d_min.
    % -----------------------------------------------------------------------
    % % Auxiliary stuff below:
    % % Plot results:
    % figure
    % subplot(1,2,1)
    % imshow(frameNoBgnd,[]);
    % hold;
    % plot(x_new,y_new,'o','Color','r','MarkerSize',5); title('obtained centre')
    % hold off;
    % % plot subarray as well:
    % subplot(1,2,2)
    % imshow(I,[]);
    % hold;
    % plot(x_new-round(x_estimate)+d+1,y_new-round(y_estimate)+d+1,'o','Color','r','MarkerSize',70); title('obtained centre')
    % hold off;
    
end