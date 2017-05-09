function [binCentres freqCounts] = fullPwD(input_vector,nbins,display_figures)
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
% Full calculation of Pair-wise Differences (PwD) for an input vector,
% calculation of their probability distribution and optional display of their
% histograms.
% 
% Note that here all pair-wise differences are calculated for all pairs, not only for
% consecutive pairs of points.
%
% Pair wise differences are calculated from start to end of input_vector:
% positive/negative pair-wise differences corresponds to values
% increasing/decreasing, respectively, from start to end of input_vector.
%
% Inputs: 
% input_vector: needs to be a COLUMN vector. 
% It is the data vector for which we want to calculate pair-wise differences.
% nbins: number of bins in histogram. (Default from Matlab is 10).
% % points_per_bin: related to the bin size for the histograms/or
% probability distribution of pair-wise differences. Default is 20 when
% there is a large unsubtracted background and around 2 when the signal has
% is background-subtracted.
% display_figures: 1 if you want to show plots, 0 if not.
%
% Outputs:
% binCentres: row vector containing the resulting bin centres of the full
% histogram of all pair-wise differences.
% freqCounts: row vector containing the occurrences (or frequency) for each
% value of pair-wise differences.
% The first plot shows the full histogram.
% The second plot shows the highest end of the full histogram, related to a
% possible large background level of the initial signal in input_vector,
% when no background subtraction has been applied.
% The third plot shows the lower end of the full histogram, related to
% photobleaching steps, blinking, noise, etc.

PD_all = allPD(input_vector); % Get column vector with all Pair-wise Differences.
% See allPD.m.

% Calculate and show pairwise difference probability distribution (plot histograms):
% nbins = round(length(PD_all)/points_per_bin); % chosen number of bins which seems to work fine in most cases. Cannot be too small!
[freqCounts,binCentres] = hist(PD_all,nbins); % histogram. 

% Display figures only if input display_figures is 1:
if display_figures == 1
    
    figure('position',[1200 100 500 1000]); % [left, bottom, width, height]; first two for lower-left corner of figure.
    subplot(3,1,1) % full histogram.
    bar(binCentres,freqCounts,'r'); % plot a bar graph of the full histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    ylim([0 1.1*max(freqCounts)]); % re-scale vertical axis.
    
    pos_negatives = find(binCentres < 0); % positions of negative bin-centre values.
    pos_positives = find(binCentres > 0); % positions of positive bin-centre values.
    if isempty(pos_negatives) % this is a fast dodgy way of sorting this out for the moment.
        pos_negatives = 1;
        warning('No negative side of histogram!');
    end
    if isempty(pos_positives) % this is a fast dodgy way of sorting this out for the moment.
        pos_positives = 1;
        warning('No positive side of histogram!');
    end
    binCentres_neg = binCentres(pos_negatives);
    freqCounts_neg = freqCounts(pos_negatives);
    binCentres_pos = binCentres(pos_positives);
    freqCounts_pos = freqCounts(pos_positives);
    
    % subplot 2:
    subplot(3,1,2) 
    bar(binCentres_neg,freqCounts_neg,'b'); % plot negative histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    title('prob distrib of I drops');
    
    % subplot 6:
    subplot(3,1,3) 
    bar(binCentres_pos,freqCounts_pos,'g'); % plot positive histogram.
    xlabel('Intensity pair-wise differences');
    ylabel('frequency');
    title(' prob distrib of I jumps > 0');
    
end