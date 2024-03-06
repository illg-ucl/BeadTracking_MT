# BeadTracking_MT - Bead Tracking on images for Magnetic Tweezer microscopy experiments

Code for the detection and tracking in 2D of spherical particles (with diffraction rings), typically used to analyse Magnetic Tweezer data. Two-dimensional (2D) tracking of beads is implemented. For tracking along the third dimension, z (beam propagation direction, optical axis), code is at the moment not included.

# Copyright and License

Copyright (c) 2017. Isabel Llorente-Garcia, Dept. of Physics and Astronomy, University College London, United Kingdom.

Released and licensed under a BSD 2-Clause License:

https://github.com/illg-ucl/BeadTracking_MT/blob/master/LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the BSD 2-Clause License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
BSD 2-Clause License for more details. You should have received 
a copy of the BSD 2-Clause License along with this program. 

Citation: If you use this software for your data analysis please acknowledge 
          it in your publications and cite as follows.
          
              -Citation example 1: 
               BeadTracking_MT software. (Version). 2017. Isabel Llorente-Garcia, 
               Dept. of Physics and Astronomy, University College London, United Kingdom.
               https://github.com/illg-ucl/BeadTracking_MT. (Download date).
               
              -Citation example 2:
               @Manual{... ,
                title  = {BeadTracking_MT software. (Version).},
                author       = {{Isabel Llorente-Garcia}},
                organization = { Dept. of Physics and Astronomy, University College London, United Kingdom.},
                address      = {Gower Place, London, UK.},
                year         = 2017,
                url          = {https://github.com/illg-ucl/BeadTracking_MT}
                }

# Guide to run code and basic steps of analysis

* Folder **"scriptsThatHelpRunCode/"** contains two Matlab scripts that help run the code. Script **"scriptToRun.m"** contains a detailed **step-by-step guide** about how to analyse a video sequence to track bead particles shown in the video. Script **"scriptToAnalyseManyVideos.m"** can be used to analyse many video files.
  
* The **basic steps** needed to track magnetic-tweezer beads on videos are as follows (more details in the above mentioned step-by-step guide):
    - Set-up Matlab and correct data folders and current directory;
    - Define "image_label" (string to automatically find the image-sequence file to analyse);
    - Set values of relevant parameters;
    - Find bead trajectories in image sequence and output them to an Excel file in the current directory using function **"FindTrajectsBeads.m"**.   
    - Link trajectory segments found into longer trajectories using **"linkTrajSegmentsBeads.m"**.
    - Each track/trajectory is identified by a Trajectory Number. Plot and save a useful .png image of the Trajectory Numbers for the found particles overlaid on top of first frame of the video (using function **"plotBeadTrajNumbers.m"**). 
    - Inspect tracks visually (on a video) and manually to validate good tracks by deciding which to accept as good using function **"goThroughBeadTracksVideo.m"**. See options for doing this quickly explained in "scriptToRun.m".
    - Analyse each track separatedly and produce one analysis Excel file and graph per track. This is based on functions **"showBeadTrajAnalysis.m"** and **"showManyBeadTrajAnalysis.m"**. The **analysis includes**: trajectory number, mean x and y positions, number of points in track, track duration, flags for good/short/long track, particle trajectory on x-y plane, first frame with trajectory overlayed, calculation of particle mean square displacement (MSD) versus delta time with error bars, fit of MSD curve (mobility and diffusion analysis), etc.
    - Note: the code does not work optimally for beads that are close together and have overlapping diffraction rings. For these, there will be poor parabolic fits to the cross-correlation peak for finding the bead centre with sub-pixel precision, so they will be excluded from
the final results later on.


# Matlab function folders

- **"openImageSequences/"**: functions to open different formats of image sequences (.sif, .dv, .tif, .m4v, .mat or .avi data) and return image data in a useful form. See **"extract_image_sequence_data.m"**. 

- **"saveDisplayOthers/"**: contains function **"saveFigurePNG.m"**, which can be used to save a current figure as a .png file.

- **"findBeadPosition/"**: functions to detect and find bead/particle positions on a single image frame. Function **"findCandidateBeadPositions.m"** finds the rough centres (x, y positions) of beads on an image frame, as candidate positions to later refine the particle-centre detection. Typically used for detecting magnetic bead positions. Magnetic beads are much larger than wavelength of illumination and their images have pronounced diffraction rings. The function recognises a dark spot at the bead centre. Function **"removeBgnd.m"** removes or substracts the background from an image of particles via two possible methods (useful for non-uniform backgrounds that would otherwise disturb subsequent steps of the analysis). Function **"excludeRegions.m"** can be used to exclude certain regions of the image. It eliminates all x-y particle positions found that fall within certain input regions (rectangular boxes). Function **"findBeadCentre1frame.m"** iteratively finds the centre of a bead on an image with sub-pixel precision, starting from a background-subtracted image (important). It extracts a cropped sub-array image of the bead, converts the bead's image to polar coordinates, calculates average radial intensity profiles along the vertical and horizontal directions for the bead, and performs cross-correlation operations to find the bead centre with sub-pixel precision (this is done for the vertical and horizontal averaged bead intensity profiles, calculating for each the cross-correlation between the bead intensity profile and its mirror image - which peaks at twice the displacement of the bead from the subarray centre - and fitting each cross-correlation function to a parabola to find the centre along x and y). Function **"eliminateCoincidentPositions.m"** eliminates near-coincident bead-centre positions (closer than a certain input distance).

- **"TrackingBeads/" - Find trajectories and ouput them to Excel file**: functions for tracking particle/bead positions on an image sequence, finding their trajectories (x, y coordinates and time) and outputting them to an Excel files. 
The two functions used are **"FindTrajectsBeads.m"** and **"linkTrajSegmentsBeads.m"**. The first one finds all particle/bead trajectories in an input image sequence, finds particle centres and joins centres in subsequent frames into trajectory segments. The second one joins segments into longer trajectories for each particle and outputs the resulting trajectory data to an Excel file. It also overlays fully connected trajectories onto the original image sequence in a video. Function **"plotBeadTrajNumbers.m"** plots Trajectory Numbers (as they appear on the Excel file generated by "linkTrajSegmentsBeads.m") next to each corresponding particle overlaid on the first frame of the chosen image sequence. This allows easy further analysis, exclusion of particles that are too close, etc.

- **"TrackingBeads/" - Inspect and validate tracks**: functions for visually inspecting and validating all particle/bead tracks found in an image sequence. Function **"goThroughBeadTracksVideo.m"** allows visual inspection and validation of the found tracks by showing a video of the found and accepted tracks (overlayed on the original bead images) so that you can label each track as good/valid (entering 1) or bad/invalid (entering 0). This process generates structures such as "good_track_nums_label.mat". This visual inspection is useful in case there are objects on the image that you want to exclude, tracks outside acceptable regions, bad tracking or other anomalies. 

- **"TrackingBeads/" - Analyse each good track**: function **"showManyBeadTrajAnalysis.m"** (which calls **"showBeadTrajAnalysis.m"**, **"getDisplacement.m"**, etc.) produces one analysis Excel file and one graph per analysed particle/bead track. A folder is generated which contains the Excel and graph files for the analysis of each track in the image sequence. It takes as input a list of "good" track numbers (previously generated). It can show the found particle trajectory overlayed on the image sequence on a video.

- **"Diffusion-MSDcalculation/"**. Contains main function **"getDisplacement.m"** which analyses trajectory data and calculates the 2D mean square displacements (MSD) as a function of time-interval delta-t for each particle trajectory (track) in a given image sequence (video).

- **"pair-wise-difference-calculation/"**: functions for the calculation of pair-wise differences (PwD) for an input vector or series of values. Used for the calculation of the Mean Square Displacement (MSD).

- Auxiliary functions. Folder "OperationsOnImages/" contains auxiliary function **"movingAvg.m"**, which can be used to calculate moving averages of two input column vectors (not really an operation on images despite the folder name). Folder "Other/" contains other auxiliary functions.
