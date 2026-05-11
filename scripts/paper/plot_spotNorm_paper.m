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
lociPath = fullfile( varPath, 'Loci SpotNorm', 'comb'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);


strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'lacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};


xNormBin = 0.04; lNormBin = 0.01; % for spotNorm plotting for all timePoints

    set( figure(1), 'Position', [600 450 340 380]) % for xNorm
    % set( figure(2), 'Position', [1010 450 400 380]) % for lNorm
    set( figure(2), 'Position', [1010 450 500 380]) % for lNorm (no folding)
    colorList = get( gca,'colororder');

c = 0;  


for j = plotNum
    
    c = c + 1;
    
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


    legtxt = sprintf( '%s (%s)%s', strain, expDate, extraName);
    legtxt2 = sprintf( '%s%s', strainName, extraName);
    legtxt = sprintf( '%s, %s', strain, strainName);


    % 1. plot abs( xNorm)
    figure(1)
    % [~, edges] = histcounts( xxNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
    %     tmp = movmean( edges, 2);   centers = tmp( 2:end);
    % errorbar( centers, mean( mX), std( mX), 'LineWidth', 2.5, 'DisplayName', legtxt, 'Color', colorList(c,:)), hold on
    [N, edges] = histcounts( xNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
        tmp = movmean( edges, 2);   centers = tmp( 2:end);
    plot( centers, N, 'LineWidth', 2.5, 'DisplayName', legtxt2), hold on


    % 2. plot lNorm
    figure(2)
    % [N, edges] = histcounts( abs( 0.5-lNorm), 'BinWidth', lNormBin, 'BinLimits', [0 0.5], 'Normalization', 'probability');
    [N, edges] = histcounts( lNorm, 'BinWidth', lNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
        tmp = movmean( edges, 2);   centers = tmp( 2:end);
    plot( centers, N, 'LineWidth', 2.5, 'DisplayName', legtxt2), hold on

            
    % display the cell numbers & spot numbers
    fprintf( '  ~~~ %3d images, %4d cells, %5d/%5d tracks,  %5d spots,  <xNorm> = %.3f    %s-%s%s\n',...
        size( cellRecord, 1), sum( cellSpots > 0), sum( cond), nTracks, numel( xNorm),...
        mean( xNorm(:), 'omitnan'), expDate, strain, extraName) % sum( cellSpots), 
    
end

%% figure setting

% xNorm
figure(1)
set( gca, 'FontName', 'Arial', 'FontSize', 20)
set( gca, 'Xtick', 0:0.2:1, 'LineWidth', 1)
% xlabel([ '|xNorm|, bin=' num2str( xNormBin)]), ylabel( 'Probability')
xlabel( '|Normalized x position|'), ylabel( 'Probability')
% title( sprintf( 'Cell with %s', extra))
legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 16)
ylim( [0 0.15])


% lNorm
figure(2)
set( gca, 'FontName', 'Arial', 'FontSize', 20)
set( gca, 'FontSize', 20, 'Xtick', 0:0.2:1, 'LineWidth', 1)
xlabel([ '|0.5-LNorm|, bin=' num2str( lNormBin)]), ylabel( 'Probability')
xlabel([ 'LNorm, bin=' num2str( lNormBin)]), ylabel( 'Probability')
xlabel( 'Normalized L position'), ylabel( 'Probability')
legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 12)
legend( 'Location', 'best', 'box', 'off', 'FontSize', 16)
% title( sprintf( 'Cell with %s', extra))
% xlim([0 0.5])
grid on


%%


outPath = 'C:\Users\yuhuanw2\Desktop\Plots';
% saveas( figure(1), fullfile(outPath, 'Fig5_xNorm.svg'));
% saveas( figure(2), fullfile(outPath, 'Fig5_LNorm.svg'));