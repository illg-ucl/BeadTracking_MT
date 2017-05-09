function frame_final = removeBgnd(frame,xx,yy,inner_radius,subarray_halfsize,method)
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
% INPUTS:
% - frame: Input image, a single frame from a recorded sequence.
% Note that the input image can be obtained as:
% frame = extract1frameB(framenumber); with framenumber=100, for example.
% Candidate bead-centre positions obtained with function findCandidateBeadPositions:
%
% - xx: columns of positions of candidate beads (horizontal). Vector.
%
% - yy: rows of positions of candidate beads (vertical). Vector.
%
% - inner_radius: radius (in pixels) of circular mask that contains each entire
% bead. Used to calculate and subtract the background: everything inside
% the circular mask is a bead, everything outside is background. Default
% ~50. Used for both methods 1 and 2.
%
% - subarray_halfsize: half-size (pixels) of subarray taken around each
% candidate bead position. It must be larger than bead radious,
% default~35-70. Used only for method 1.
%
% - method: enter 1 or 2 as below (default = 1):
% Method 1 (fast one): Take the image frame, exclude pixels within the
% circular masks (beads masks) at the detected positions (xx,yy), then take
% a square subarray around the circular masks and calculate the average
% local background intensity per pixel in the background-only pixels (outside
% circular masks) within the subarrays. For the subarray around each bead,
% subtract the average background intensity from all pixels around that
% bead in the original frame.
% Method 2 (slow one): Take the image frame, exclude pixels within the
% circular masks (beads masks) at the detected positions (xx,yy), then
% interpolate to obtain background within beads regions and obtain an
% all-background frame. Then subtract the all-background frame from the
% original frame. It takes a few seconds to run for a ~1024x1024 image on a
% decent computer!
%
% OUTPUTS:
% - frame_final: frame with background subtracted.
%
% Example of how to run this function: frame1_final = removeBgnd(frame1,x1,y1,50,60,1);

% Preallocate a matrix of zeros:
bead_mask = zeros(size(frame,1),size(frame,2));

% Create matrices containing the x and y positions corresponding to the
% frame. Xs is a matrix of the same size as frame, containing y values for
% all pixels and similarly for Ys. In Xs, x varies along second dimension
% (horizontal) in matrix. In Ys, y varies along first dimension (vertical).
[Xs,Ys] = meshgrid(1:size(frame,2),1:size(frame,1));

% Loop through all candidate bead positions:
for i = 1:length(xx)
    % First assign value 1 to pixels where beads are, within circular masks, for
    % which radius <= inner_radius. C = hypot(A,B) returns SQRT(ABS(A).^2+ABS(B).^2).
    bead_mask = bead_mask | hypot(Xs-xx(i), Ys-yy(i)) <= inner_radius;   % Beads mask (binary image).
end

% Background mask is the negative of the beads mask:
bgnd_mask = double(~bead_mask); % background mask (binary image).
Ibgnd = frame.*bgnd_mask; % bgnd-only image with holes at bead positions.

%% Method 1 (fast): subtract local average background intensity per pixel:
if method == 1
    % Preallocate final frames:
    frame_bgnd = Ibgnd; % the holes at bead positions are later filled in with average background values.
    % Loop through all candidate bead positions:
    for i = 1:length(xx)       
        % For selecting subarray around each candidate bead position, initally
        % the half-sizes of the subarray along each direction are equal to the
        % input subarray_halfsize. They are only reassigned if a bead
        % is at the edge of the image.
        d_top = subarray_halfsize;
        d_bottom = subarray_halfsize;
        d_left = subarray_halfsize;
        d_right = subarray_halfsize;
        % If the bead is at edge of image, take a smaller subarray around
        % it but as large as possible until the edge is reached,
        % so reasign the d values:
        if (round(yy(i)) - d_top) < 1
            d_top = round(yy(i)) - 1;
        end
        if (round(yy(i)) + d_bottom) > size(frame,1)
            d_bottom = size(frame,1) - round(yy(i));
        end
        if (round(xx(i)) - d_left) < 1
            d_left = round(xx(i)) - 1;
        end
        if (round(xx(i)) + d_right) > size(frame,2)
            d_right = size(frame,2) - round(xx(i));
        end
        % Create image subarray centered on bead position:
        frame_sub = frame(round(yy(i))-d_top:round(yy(i))+d_bottom,round(xx(i))-d_left:round(xx(i))+d_right);
        % Create background mask subarray centered on bead position (binary subarray):
        bgnd_mask_sub = bgnd_mask(round(yy(i))-d_top:round(yy(i))+d_bottom,round(xx(i))-d_left:round(xx(i))+d_right);
        % Calculate average background per pixel around bead position in
        % background area:
        pos_bgnd_sub = find(bgnd_mask_sub==1); % vector of positions of all background pixels in background mask subarray.
        I_bgnd_vector_sub = frame_sub(pos_bgnd_sub); % vector of background intensity values.
        Ibgnd_avg = mean(I_bgnd_vector_sub); % mean background value per pixel around bead.
        % Asign the mean average background intensity to all pixels within circular area occupied by bead:   
        frame_bgnd(hypot(Xs-xx(i),Ys-yy(i)) <= inner_radius) = Ibgnd_avg; % Background-only (everywhere) image.
    end
    % Final background-subtracted frame:
    frame_final = frame - frame_bgnd;
    
    % Aid plots:
%     figure;
%     subplot(2,2,1)
%     imshow(frame,[]);
%     title('Original frame')
%     subplot(2,2,2)
%     imshow(Ibgnd,[])
%     title('Background region')
%     subplot(2,2,3)
%     imshow(frame_bgnd,[])
%     title('Background-only frame')
%     subplot(2,2,4)
%     imshow(frame_final,[])
%     title('Background-subtracted frame')
end


%% Method 2 (slow): interpolate background where beads are:
if method == 2
    
    pos_beads = find(bead_mask==1); % positions of beads.
    x_beads = Xs(pos_beads); % column vector of x-positions of bead regions.
    y_beads = Ys(pos_beads); % column vector of y-positions of bead regions.
    
    pos_bgnd = find(bgnd_mask==1); % positions of all background pixels.
    x_bgnd = Xs(pos_bgnd); % vector of x-positions of bgnd regions.
    y_bgnd = Ys(pos_bgnd); % vector of y-positions of bgnd regions.
    I_bgnd_vector = Ibgnd(pos_bgnd); % vector of background intensity values.
    
    % Interpolate to find background at positions of beads:
    Ibgnd_interpolated = griddata(x_bgnd,y_bgnd,I_bgnd_vector,x_beads,y_beads,'natural');
    % Method: 'natural' is the only method that does not introduce lines.
    
    % Complete background image (as a vector):
    x_fullBgnd = [x_beads; x_bgnd];
    y_fullBgnd = [y_beads; y_bgnd];
    I_fullBgnd = [Ibgnd_interpolated; I_bgnd_vector];
    
    % Calculate background-corrected subarray image:
    frame_bgnd = zeros(size(frame,1),size(frame,2)); % bgnd-only frame
    frame_final = zeros(size(frame,1),size(frame,2)); % bgnd-subtracted frame
    for i = 1:length(x_fullBgnd)
        x = x_fullBgnd(i);
        y = y_fullBgnd(i);
        bgnd_intensity = I_fullBgnd(i);
        frame_bgnd(y,x) = bgnd_intensity;
        frame_final(y,x) = frame(y,x)-bgnd_intensity;
    end
    
%     % Aid plots:
%     figure;
%     subplot(2,3,1)
%     imshow(frame,[]);
%     title('Original frame')
%     subplot(2,3,2)
%     imshow(Ibgnd,[])
%     title('Background region')
%     subplot(2,3,3)
%     plot(x_beads,-y_beads,'.');
%     xlim([1 size(frame,2)])
%     ylim([-size(frame,1) 1])
%     axis('equal')
%     title('Bead positions')
%     subplot(2,3,4)
%     plot(x_bgnd,-y_bgnd,'.');
%     xlim([1 size(frame,2)])
%     ylim([-size(frame,1) 1])
%     axis('equal')
%     title('Background positions')
%     subplot(2,3,5)
%     imshow(frame_bgnd,[])
%     title('Background-only frame')
%     subplot(2,3,6)
%     imshow(frame_final,[])
%     title('Background-subtracted frame')
    
end

% % Aid plots:
% figure;
% subplot(1,2,1)
% imshow(frame_bgnd,[])
% subplot(1,2,2)
% imshow(frame_final,[])

end