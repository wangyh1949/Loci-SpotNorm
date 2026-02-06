%{
Author: Yu-Huan Wang   (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 2/6/2026
    Last update date: 2/6/2026

    ~~~~~~ adapted from sim_xNorm_loci.m ~~~~~~

Description: This code numerically simulate the xNorm of loci with a DNA
free zone. 

It plots xNorm of loci with different width with a DNA free zone whose
radius is r.

I assumed loci position is homogeneously distributed in a cylinder
(effectively like cytoplasmic proteins)
---------------------------------------------------------------------------
%}


clear, clc, close all
 
locErr = 20; 

% scaling factor for locErr
scaleF = 2; % here assumes cell has R = 1 um, which is ~2x typical cell size

lociWidList = [0.3 0.35 0.4 0.45 0.5 0.6 0.7]; % width of loci position

figure( 'Position', [600 500 400 380]), hold on
colorList = flip( winter( numel( lociWidList)));


c = 0;

for lociWid = lociWidList

    c = c + 1;

    legtxt = sprintf( 'wid = %.2f μm', lociWid);

    R = lociWid; % radius of loci position

    % r = 0; % radius of DNA excluded region
    r = R - 0.2;

    % numerical xNorm setup
    nBin = 100;    binEdges = linspace( -R, R, nBin+1);
    tmp = movmean( binEdges, 2);    binCenters = tmp(2:end);
    
    % [Membrane]     calculate the surface density (arc length) in each bin
    surf = -diff( acos( binEdges)); % arccos( x/r) = theta, r* theta = arc, r=1
    
    % [Cytoplasmic]  calculate the area in each bin
    cyto = sqrt( R^2- binCenters.^2); % corresponding y of each bin centers

        % calculate excluded region
        cyto_exclude = cyto - sqrt( r^2 - binCenters.^2);
        region_exclude = abs( binCenters) < r;
        cyto( region_exclude) = cyto_exclude( region_exclude); 

    % Gaussian blue vector [1,51] in x axis, the locErr was scaled by scaleF
    blurX = normpdf( binEdges, 0, locErr* scaleF/ 1000);

    % apply locErr in x axis
    surfLoc = conv( surf, blurX);
    cytoLoc = conv( cyto, blurX);

    % conv: [1,50]*[1,51]-->[1,100], new binEdges:[-2,2] with 2*nBin
    newX = linspace( binCenters(1)-R, binCenters(end)+R, nBin*2); % x coordinate after convolution

    % interpolate for the same coordinate
    xNum = 0: 0.01: 1;
    vq = interpn( newX, cytoLoc, xNum, 'cubic'); % interpolated value at the same bin as experimental xNorm
    normF = sum( vq, 'omitnan');    xNormNum = vq/ normF* 4; % numerical calculated xNorm
    xNormNum( abs( xNormNum) < 1e-10) = 0; % get rid of tiny negative value from interp

    % plot xNorm
    plot( xNum, xNormNum, 'LineWidth', 2.5, 'color', colorList( c,:), ...
        'DisplayName', legtxt)
end

set( gca, 'FontSize', 14, 'Xtick', 0:0.2:1, 'LineWidth', 1)
legend( 'location', 'northeast', 'FontSize', 11, 'box', 'off')
xlabel( 'xNorm'), ylabel( 'Probability'), box on
title( 'Numerical xNorm', 'FontSize', 16)
ylim( [0 0.15])