%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 12/22/2025
    Last update date: 12/22/2025

~~~~~~~~~~~ adapted from plotSpotNorm_loci.m ~~~~~~~~~~~
    for loci imaging experiments (loci position)

Description: This script plot xNorm & lNorm for loci tracking data
---------------------------------------------------------------------------
%}


clear, clc, close all

% set the path to store analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);


strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'LacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};


xNormBin = 0.04; lNormBin = 0.02; % for spotNorm plotting for all timePoints
xNormStat = nan( length( plotNum), 2);
labelName = cell( length( plotNum), 1);

    set( figure(1), 'Position', [600 450 400 380]) % for xNorm
    set( figure(2), 'Position', [1010 450 400 380]) % for lNorm
    colorList = get( gca,'colororder');

c = 0;  


for j = plotNum
    
    c = c + 1;
    
    % load lociPos file
    load( fullfile( lociList( j).folder, lociList( j).name))

    index = find( strcmp( string( strainList), strain(3:end)));
    if isempty( index), strainName = [ strain extraName];
    else, strainName = [ nameList{ index} extraName]; end


    % count spot number for each cell
    cellSpots = accumarray( cellNum(:), 1, [totalCells, 1]);

    % ~~~~~~~~~ condition for cell variables ~~~~~~~~~
    goodCells = find( cellSpots > 0);   extra = '1+ spot'; 
    % goodCells = find( cellSpots == 1);   extra = '1 spot';
    % goodCells = find( cellSpots == 2);   extra = '2 spots';


    xxNorm = tracksxNorm(:,1);   xNorm40 = tracksxNorm40;
    llNorm = tracksLNorm(:,1);   lNorm40 = tracksLNorm40;
    mid = tracksMid(:,1);       mid40 = tracksMid40;
    % tracksMid: [first spot, whole track]


    % ~~~~~~~~ Condition ~~~~~~~~~
    condSpots = ismember( cellNum, goodCells); % flag for spots in selected cells
    condxNorm = min( mid40, [], 2); % exclude tracks with any cap points
    cond = condSpots & condxNorm;

    % xNorm = abs( xxNorm( cond));
    % lNorm = llNorm( cond);
    xNorm = abs( xNorm40( cond,:));
    lNorm = lNorm40( cond,:);


    legtxt = sprintf( '%s (%s)', strain, expDate);
    legtxt2 = sprintf( '%s%s', strainName);


    % 1. plot abs( xNorm)
    figure(1)
    % [~, edges] = histcounts( xxNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
    %     tmp = movmean( edges, 2);   centers = tmp( 2:end);
    % errorbar( centers, mean( mX), std( mX), 'LineWidth', 2.5, 'DisplayName', legtxt, 'Color', colorList(c,:)), hold on
    [N, edges] = histcounts( xNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
        tmp = movmean( edges, 2);   centers = tmp( 2:end);
    plot( centers, N, 'LineWidth', 2.5, 'DisplayName', legtxt), hold on


    % 2. plot lNorm
    figure(2)
    [N, edges] = histcounts( abs( 0.5-lNorm), 'BinWidth', lNormBin, 'BinLimits', [0 0.5], 'Normalization', 'probability');
        tmp = movmean( edges, 2);   centers = tmp( 2:end);
    plot( centers, N, 'LineWidth', 2.5, 'DisplayName', legtxt2), hold on

            
    % display the cell numbers & spot numbers
    fprintf( '  ~~~ %3d images,  %4d cells, %5d/%5d spots,  <xNorm> = %.3f,    %s-%s%s\n',...
        size( cellRecord, 1), sum( cellSpots > 0), length( xNorm),...
        sum( cellSpots), mean( xNorm, 'omitnan'), expDate, strain, extraName) % sum( cellSpots), 
    
end

%% figure setting

% xNorm
figure(1)
set( gca, 'FontSize', 14, 'Xtick', 0:0.2:1, 'LineWidth', 1)
xlabel([ '|xNorm|, bin=' num2str( xNormBin)]), ylabel( 'Probability')
title( sprintf( 'Cell with %s', extra))
legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 11)
ylim( [0 0.13])

% lNorm
figure(2)
set( gca, 'FontSize', 14, 'Xtick', 0:0.2:1, 'LineWidth', 1)
xlabel([ '|0.5-LNorm|, bin=' num2str( lNormBin)]), ylabel( 'Probability')

legend( 'Location', 'southeast', 'box', 'off', 'FontSize', 11)
title( sprintf( 'Cell with %s', extra))
xlim([0 0.5])
