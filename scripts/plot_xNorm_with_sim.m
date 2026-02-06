%{
Author: Yu-Huan Wang   (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 2/6/2026
    Last update date: 2/6/2026

    ~~~~~~ adapted from sim_xNorm_loci.m & plot_spotNorm.m ~~~~~~

Description: This code plot the overlay of loci xNorm with the numerically
simulated xNorm

I assumed loci position is homogeneously distributed in a cylinder
(effectively like cytoplasmic proteins), 
---------------------------------------------------------------------------
%}


clear, clc, close all

% parameters
locErr = 20; wid = 1;

lociWidList = [0.35 0.4 0.45 0.5 0.55 0.6]; % width of loci

% figure( 'Position', [600 500 400 380]), hold on

% numerical xNorm setup
    nBin = 100;    binEdges = linspace( -1, 1, nBin+1);
    tmp = movmean( binEdges, 2);    binCenters = tmp(2:end);    
    surf = -diff( acos( binEdges)); % membrane
    y = sqrt( 1- binCenters.^2); % cytoplasmic


%% Plot experimental xNorm

% set the path to store analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
% lociPath = fullfile( varPath, 'Loci SpotNorm', 'single day'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);


strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'LacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};

xNormBin = 0.04; % for spotNorm plotting for all timePoints

cc = 0;  

for j = plotNum
    
    cc = cc + 1;
    

    % 1. plot simulated xNorm

        figure( 'Position', [600+30*cc 500+30*cc 400 380]), hold on
        c = 0;  colorList = flip( winter( numel( lociWidList)));
    
        for lociWid = lociWidList
        
            c = c + 1;        
            legtxt = sprintf( 'wid = %.2f μm', lociWid);    
        
            dilF = wid/ lociWid; % calculate dilF
            R = wid/2;  scaleF = dilF/ R; % rescale model cell
            blurX = normpdf( binEdges, 0, locErr* scaleF/ 1000); % Gaussian blur vector
        
            % apply locErr in x axis
            cytoLoc = conv( y, blurX);
        
            % conv: [1,50]*[1,51]-->[1,100], new binEdges:[-2,2] with 2*nBin
            conX = linspace( binCenters(1)-1, binCenters(end)+1, nBin*2); % x coordinate after convolution
            newX = conX/ dilF; % real bin centers after scaling by dilF (membrane/outline)
        
            % interpolate for the same coordinate
            xNum = 0: 0.01: 1;
            vq = interpn( newX, cytoLoc, xNum, 'cubic'); % interpolated value at the same bin as experimental xNorm
            normF = sum( vq, 'omitnan');    xNormNum = vq/ normF* 4; % numerical calculated xNorm
        
            % plot xNorm
            plot( xNum, xNormNum, 'LineWidth', 2, 'color', colorList( c,:), ...
                'DisplayName', legtxt)
        end


    % load lociPos file
    load( fullfile( lociList( j).folder, lociList( j).name))

    % find strain name
    index = find( strcmp( string( strainList), strain(3:end)));    
    if isempty( index), strainName = strain;
    else, strainName = nameList{ index}; end

        tInt = findInt( extraName, expDate, strain);
        extraName = erase( extraName, [ ' ' tInt]);
        expDate = erase( expDate, ' comb');


    % count spot number for each cell
    cellSpots = accumarray( cellNum(:), 1, [totalCells, 1]);

    % ~~~~~~~~~ condition for cell variables ~~~~~~~~~
    goodCells = find( cellSpots > 0);   extra = '1+ spot'; 
    % goodCells = find( cellSpots == 1);   extra = '1 spot';
    % goodCells = find( cellSpots == 2);   extra = '2 spots';


    xxNorm = tracksxNorm(:,1);   xNorm40 = tracksxNorm40;
    mid = tracksMid(:,1);       mid40 = tracksMid40;
    % tracksMid: [first spot, whole track]


    % ~~~~~~~~ Condition ~~~~~~~~~
    condSpots = ismember( cellNum, goodCells); % flag for spots in selected cells
    condxNorm = min( mid40, [], 2); % exclude tracks with any cap points
    cond = condSpots & condxNorm;

    % xNorm = abs( xxNorm( cond));
    xNorm = abs( xNorm40( cond,:));

    legtxt = sprintf( '%s (%s)%s', strain, expDate, extraName);
    legtxt2 = sprintf( '%s%s', strainName, extraName);
    legtxt = sprintf( '%s %s', strain, strainName);
    

    % 2. plot abs( xNorm)
    [N, edges] = histcounts( xNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
        tmp = movmean( edges, 2);   centers = tmp( 2:end);
    plot( centers, N, 'r', 'LineWidth', 2.5, 'DisplayName', legtxt, 'HandleVisibility', 'off'), hold on

    % figure setting
    set( gca, 'FontSize', 14, 'Xtick', 0:0.2:1, 'LineWidth', 1)
    xlabel( sprintf( 'xNorm (locE=%dnm)', locErr)), ylabel( 'Probability'), box on
    legend( 'location', 'northeast', 'box', 'off', 'FontSize', 11)
    title( legtxt)
    ylim( [0 0.15])
end