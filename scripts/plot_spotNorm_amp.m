%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 3/27/2026
    Last update date: 3/27/2026

~~~~~~~~~~~ adapted from plot_spotNorm.m & plot_subMSD_amp.m ~~~~~~~~~~~


Description: This script plot xNorm & lNorm for loci tracking data binned
by signal amplitude
---------------------------------------------------------------------------
%}


clear, clc, close all

% set the path to store analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
% lociPath = fullfile( varPath, 'Loci SpotNorm', 'comb'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);


strainList = { 727 728 729 730 731 734 701 662 311 725 830};
nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' '12tetO@LacZ' ... 'LacZ'...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ' 'sal'};


xNormBin = 0.04; lNormBin = 0.01; % for spotNorm plotting for all timePoints


for j = plotNum
    
    % set up figures for each strain
    f1 = figure; set( f1, 'Position', [600 450 400 380]) % for xNorm
    f2 = figure; set( f2, 'Position', [1010 450 400 380]) % for lNorm
    
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


    % ~~~~~~~~ Condition ~~~~~~~~~
    mid = tracksMid( :, 1);        mid40 = tracksMid40;
    condxNorm = min( mid40, [], 2); % exclude tracks with any cap points
    condSpots = ismember( cellNum, goodCells); % flag for spots in selected cells

    cond = condSpots & condxNorm; % flag for tracks with selected spots and good xNorm

    xxNorm = tracksxNorm( cond, 1);   xNorm40 = tracksxNorm40( cond,:);
    llNorm = tracksLNorm( cond, 1);   lNorm40 = tracksLNorm40( cond,:);
    % tracksMid: [first spot, whole track]
    

    % divide data by signal amplitude
    binData = tracksAmp( cond, 1);  binName = 'amp'; % use amplitude at first frame for binning
    % binData = mean(tracksAmp(:,1:40), 2);  binName = 'amp40'; % use average amplitude of first 40 frames for binning

    binPer = [0 0.1 0.3 0.5 0.7 0.9 1]; % binning percentage
    binList = quantile( binData, binPer); 
    colorList = flip( winter( numel( binPer))); c = 0; 

                
        % display the cell numbers & spot numbers
        fprintf( '  ~~~ %3d images, %4d cells, %5d/%5d tracks,  %5d spots,  <xNorm> = %.3f    %s-%s%s\n',...
            size( cellRecord, 1), sum( cellSpots > 0), sum( cond), nTracks, numel( xxNorm),...
            mean( xxNorm(:), 'omitnan'), expDate, strain, extraName)

    for k = 1: numel( binList)- 1

        c = c + 1;
        % ~~~~~~~~ Condition ~~~~~~~~~
        condAmp = binData >= binList(k) & binData < binList(k+1); % amp in this bin

        % xNorm = abs( xxNorm( cond));
        % lNorm = llNorm( cond);
        xNorm = abs( xNorm40( condAmp,:));
        lNorm = lNorm40( condAmp,:);

        % legtxt = sprintf( '%s: %g%%-%g%%', binName, binPer(k)*100, binPer(k+1)*100);
        legtxt = sprintf( '%s: %g%%', binName, binPer(k+1)*100);
        fprintf( '      %d/%d tracks,   %s\n', sum( condAmp), sum( cond), legtxt)


        % 1. plot abs( xNorm)
        figure(f1), hold on
        % [~, edges] = histcounts( xxNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
        %     tmp = movmean( edges, 2);   centers = tmp( 2:end);
        % errorbar( centers, mean( mX), std( mX), 'LineWidth', 2.5, 'DisplayName', legtxt, 'Color', colorList(c,:)), hold on
        [N, edges] = histcounts( xNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
            tmp = movmean( edges, 2);   centers = tmp( 2:end);
        plot( centers, N, 'LineWidth', 2.5, 'Color', colorList(c,:), 'DisplayName', legtxt)


        % 2. plot lNorm
        figure(f2), hold on
        [N, edges] = histcounts( abs( 0.5-lNorm), 'BinWidth', lNormBin, 'BinLimits', [0 0.5], 'Normalization', 'probability');
        % [N, edges] = histcounts( lNorm, 'BinWidth', lNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
            tmp = movmean( edges, 2);   centers = tmp( 2:end);
        plot( centers, N, 'LineWidth', 2.5, 'Color', colorList(c,:), 'DisplayName', legtxt)

    end

    
    % figure setting
    titletxt = sprintf( '%s, %s%s', strain, strainName, extraName);

    % xNorm
    figure(f1)
    set( gca, 'FontSize', 14, 'Xtick', 0:0.2:1, 'LineWidth', 1)
    xlabel([ '|xNorm|, bin=' num2str( xNormBin)]), ylabel( 'Probability')
    title( titletxt)
    legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 11)
    ylim( [0 0.15])

    % lNorm
    figure(f2)
    set( gca, 'FontSize', 14, 'Xtick', 0:0.2:1, 'LineWidth', 1)
    xlabel([ '|0.5-LNorm|, bin=' num2str( lNormBin)]), ylabel( 'Probability')
    % xlabel([ 'LNorm, bin=' num2str( lNormBin)]), ylabel( 'Probability')
    legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 11)
    % legend( 'Location', 'best', 'box', 'off', 'FontSize', 12)
    title( titletxt)
    xlim([0 0.5])
    % grid on
end
