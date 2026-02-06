%{
Author: Yu-Huan Wang   (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 2/3/2026
    Last update date: 2/3/2026

    ~~~~~~ adapted from xNormSimPlotTest.m ~~~~~~

Description: This code numerically simulate the xNorm of loci

I assumed loci position is homogeneously distributed in a cylinder
(effectively like cytoplasmic proteins)
---------------------------------------------------------------------------
%}


clear, clc, close all

% parameters
locErr = 20; wid = 1;

lociWidList = [0.35 0.4 0.45 0.5 0.55]; % width of loci

figure( 'Position', [600 500 400 380]), hold on
colorList = flip( winter( numel( lociWidList)));

% numerical xNorm setup
    nBin = 100;    binEdges = linspace( -1, 1, nBin+1);
    tmp = movmean( binEdges, 2);    binCenters = tmp(2:end);
    
    % [Membrane]     calculate the surface density (arc length) in each bin
    surf = -diff( acos( binEdges)); % arccos( x/r) = theta, r* theta = arc, r=1
    
    % [Cytoplasmic]  calculate the area in each bin
    y = sqrt( 1- binCenters.^2); % corresponding y of each bin centers

c = 0;

for lociWid = lociWidList

    c = c + 1;

    legtxt = sprintf( 'wid = %.2f μm', lociWid);

    % calculate dilF (effective ratio of cell width & loci width)
    dilF = wid/ lociWid;

    % rescale the model cell to match real inner membrane,  r=1 --> R/dilF
    R = wid/2;  scaleF = dilF/ R;

    % Gaussian blue vector [1,51] in x axis, the locErr was scaled by scaleF
    blurX = normpdf( binEdges, 0, locErr* scaleF/ 1000);

    % apply locErr in x axis
    surfLoc = conv( surf, blurX);
    cytoLoc = conv( y, blurX);

    % conv: [1,50]*[1,51]-->[1,100], new binEdges:[-2,2] with 2*nBin
    conX = linspace( binCenters(1)-1, binCenters(end)+1, nBin*2); % x coordinate after convolution

    newX = conX/ dilF; % real bin centers after scaling by dilF (membrane/outline)

    % % plot xNorm
    % plot( newX, cytoLoc/ sum( cytoLoc), 'LineWidth', 3, 'color', colorList( c,:), ...
    %     'DisplayName', legtxt)
    % plot( newX, surfLoc/ sum( surfLoc), 'LineWidth', 3, 'color', colorList( c,:), ...
    %     'DisplayName', legtxt)

    % interpolate for the same coordinate
    xNum = 0: 0.01: 1;
    vq = interpn( newX, cytoLoc, xNum, 'cubic'); % interpolated value at the same bin as experimental xNorm
    normF = sum( vq, 'omitnan');    xNormNum = vq/ normF* 4; % numerical calculated xNorm

    % plot xNorm
    plot( xNum, xNormNum, 'LineWidth', 2.5, 'color', colorList( c,:), ...
        'DisplayName', legtxt)
end

set( gca, 'FontSize', 14, 'Xtick', 0:0.2:1, 'LineWidth', 1)
legend( 'location', 'northeast', 'FontSize', 12, 'box', 'off')
xlabel( 'xNorm'), ylabel( 'Probability'), box on
title( 'Numerical xNorm', 'FontSize', 16)
ylim( [0 0.15])