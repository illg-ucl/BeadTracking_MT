function [xx,yy] = findCandidateBeadPositions(frame,method)
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
% Find candidate x, y positions of magnetic beads on images. Magnetic beads
% are much larger than wavelength of illumination and their images have
% pronounced diffraction rings. The function recognises dark spot in
% centre. Two methods are used: (1) morphological operations to obtain
% candidates for bead centres; (2) Centre of Mass method.
% Note that the first method appears to be more accurate. The second method
% is far more sensitive to uneven background illumination and nearby beads.
% None of this methods give sub-pixel precision, though.
%
% The candidate spots can later be used as input of function
% "findBeadCentre1frame".
%
% INPUTS:
% - frame: Input image, a single frame from a recorded sequence. 
% Note that the input image can be obtained as:   
% frame = extract1frameB(framenumber); with framenumber=100, for example.
% - method: 1 for  morphological operations to obtain candidates for bead 
%   centres (thresholding, eroding, etc.); 
%           2 for centre-of-mass method;
%           0 to plot and compare both methods (output is that of Method 1);
%
% OUTPUTS:
% Candidate bead-centre positions obtained from chosen method:
% - xx: columns of positions of spot candidates (horizontal).
% - yy: rows of positions of spot candidates (vertical).
%
%
% It can be useful to apply a Gaussian convolution to filter the image before
% using this function (e.g.  frame2 = imfilter(frame,fspecial('gaussian',3))).


%% PARAMETERS:
disk_radius1 = 3; % disk radius in pixels, for smoothing transformations below. Default 4.
disk_radius2 = 1; % disk radius in pixels for eliminating small things that are not beads. Default 6.
threshold_factor = 0.72; % For thresholding, aiming to detect the black centre of magnetic beads. Default 0.7.
% Adjust according to image contrast.
% A factor ~0.6 seems to be best to allow detection of dimmer beads and also separation
% of nearby beads, though it can give rise to false candidates (between close beads).
% A factor ~0.4 times the automatic threshold helps separate nearby beads
% but does not allow the detection of dimmer beads. 
aroundBead = 35; % For centre of mass method: halfsize of subarray (in pixels) taken around each bead. Default 35-50.


%% Method 1: morphological operations:

% Rescale input image frame to have values between 0 and 1: 
frame0 = mat2gray(frame); 

se1 = strel('disk',disk_radius1); % structural element, disk of radius disk_radius1 (see PARAMETERS).
se2 = strel('disk',disk_radius2); % structural element, disk of radius disk_radius2 (see PARAMETERS).

% Smooth out the background with a mrphological opening operation (= erosion
% followed by dilation) with se:
frameSmooth = imopen(frame0,se1);

% Thresholding: 
% graythresh gives the threshold and im2bw converts image into a
% thresholded image with only black (zeros) and white (ones) pixels.
% BgndMask: background mask: 0 at bead regions, 1 in rest.
% SignalMask: 1 at bead regions, 0 in rest.
BgndMask = im2bw(frameSmooth,threshold_factor*graythresh(frameSmooth));
SignalMask = imcomplement(BgndMask);

% This eliminates small pixel areas which are not actually beads:
SignalMask2 = imopen(SignalMask,se2); 

% Reduce beads to single points:
% "ultimate erosion" of a binary image = regional maxima of the Euclidean
% distance transform of the complement of the image.
SignalMask3 = bwulterode(SignalMask2);
% Shrink objects to single points:
SignalMask4 = bwmorph(SignalMask3,'shrink',1);
% In principle, sum(sum(signalMask4)) is the number of beads in the image.


% Create matrices containing the x and y positions corresponding to the
% subarray I. Xpos is a matrix of the same size as I, containing x values for
% all pixels and similarly for Ypos.
[Xpos Ypos] = meshgrid(1:size(frame,2),1:size(frame,1));

% The positions of the candidate bead positions are:
xx0 = nonzeros(Xpos.*SignalMask4); % "nonzeros" returns a column vector with the non-zero matrix elements ordered by columns.
yy0 = nonzeros(Ypos.*SignalMask4);


if method ==1
    % Generate plot of results for Method 1 (morphological operations above):
%     figure
%     imshow(frame0,[],'InitialMagnification',150);
%     hold;
%     % plot candidates from morphological method as red circles:
%     plot(xx0,yy0,'o','Color','r','MarkerSize',12); title('obtained centres')
%     hold off;
    
    % Generate output:
    xx = xx0;
    yy = yy0;
end


%% Method 2. Centre-of-mass method
% Method 2 needs method 1.

if method == 2 || method == 0
    % Use CENTRE OF MASS METHOD to obtain slightly more precise bead centres
    % (though still not with sub-pixel precision).
    
    % Halfsize of subarray (in pixels) taken around each bead candidate
    % position (see PARAMETERS above):
    d = aroundBead;
    
    % Preallocate new vectors of refined positions:
    xx_com = zeros(length(xx0),1);
    yy_com = zeros(length(yy0),1);
    
    % Loop through all candidate positions found:
    for i = 1:length(xx0)
        
        xx_0 = xx0(i);
        yy_0 = yy0(i);
        % If the spot candidate is at the edge of the image move away from the edge
        % until we can take a subarray around it:
        if (round(yy_0)-d)<1
            yy_0 = d+1;
        end
        if (round(yy_0)+d)>size(frame,1)
            yy_0 = size(frame,1)-d;
        end
        if (round(xx_0)-d)<1
            xx_0 = d+1;
        end
        if (round(xx_0)+d)>size(frame,2)
            xx_0 = size(frame,2)-d;
        end
        
        % Create image subarray (Isub) of size (2*d+1)x(2*d+1) centred on the
        % centroid estimate pixel (xx(i),yy(i)). I is fixed during the iterative process of finding the centre of the spot.
        Isub = frame0(round(yy_0)-d:round(yy_0)+d,round(xx_0)-d:round(xx_0)+d);
        
        Isub2 = Isub;
        % Subtract mean and take absolute value - actually, this works worse if
        % there are gradients in the image background, etc, so omit.
        % Isub2 = abs(Isub-mean2(Isub)); % re-scaled intensity to avoid
        % negative values for centre-of-mass method.
        
        % Create matrices containing the x and y positions corresponding to the
        % subarray Isub. Xs is a matrix of the same size as Isub, containing x values for
        % all pixels and similarly for Ys.
        [Xs,Ys] = meshgrid(round(xx_0)-d:round(xx_0)+d,round(yy_0)-d:round(yy_0)+d);
        
        % Calculate revised estimates of x and y bead centre position by centre
        % of mass method: weighing positions by the re-scaled intensity Isub2:
        xx_com(i) = sum(sum(Isub2.*Xs))/sum(sum(Isub2));
        yy_com(i) = sum(sum(Isub2.*Ys))/sum(sum(Isub2));
        
    end
    
%     % Generate plot of results:
%     figure
%     imshow(frame0,[],'InitialMagnification',150);
%     hold;
%     % plot candidates from com method as green circles:
%     plot(xx_com,yy_com,'o','Color','g','MarkerSize',14);
%     hold off;
    
    % Generate output:
    xx = xx_com;
    yy = yy_com;
    
end

if method == 0
    % Compare results of Method 1 and Method 2:
    figure
    imshow(frame0,[],'InitialMagnification',150);
    hold;
    % plot candidates from morphological method as red circles:
    plot(xx0,yy0,'o','Color','r','MarkerSize',12); title('obtained centres')
    % plot candidates from com method as green circles:
    plot(xx_com,yy_com,'o','Color','g','MarkerSize',14);
    hold off;
    
    % Generate output (results of Method 1 (arbitrary)):
    xx = xx0;
    yy = yy0;
end


end

