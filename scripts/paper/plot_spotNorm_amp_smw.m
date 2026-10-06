%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 7/14/2026
    Last update date: 7/14/2026

~~~~~~~~~~~ adapted from plot_spotNorm_amp.m ~~~~~~~~~~~

Description: This script plot xNorm & lNorm for loci tracking data binned
by signal amplitude

this script is specifically for SK830
    Loci oufti SK830 260520 comb 20-200ms.mat
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

cc = 0;

j = plotNum;
    
cc = cc + 1;
% set up figures for each strain
f1 = figure( 'Position', [500 400 400 370]); % for xNorm
f3 = figure( 'Position', [910 450 400 280]);

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
if ~exist( 'tracksSig', 'var')
    sigFactor = 1; % if no signal info, just use amplitude
else
    sigFactor = tracksSig( cond); % use signal at first frame for binning
end
binData = tracksAmp( cond,1).* sigFactor;  binName = 'amp'; % use amplitude at first frame for binning
% binData = mean(tracksAmp(:,1:40), 2);  binName = 'amp40'; % use average amplitude of first 40 frames for binning

% binPer = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; % SK830
binPer = [0 0.02 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]; % SK830
binList = quantile( binData, binPer); 
colorList = flip( winter( numel( binPer))); c = 0; 
            
    % display the cell numbers & spot numbers
    fprintf( '  ~~~ %3d images, %4d cells, %5d/%5d tracks,  %5d spots,  <xNorm> = %.3f    %s-%s%s\n',...
        size( cellRecord, 1), sum( cellSpots > 0), sum( cond), nTracks, numel( xxNorm),...
        mean( xxNorm(:), 'omitnan'), expDate, strain, extraName)


% 1. plot intensity distribution
figure( f3)
histogram( log( binData), 50, ... 'binWidth', 0.01,
    'EdgeColor', 'none', 'FaceAlpha', 0.5)
xline( log( binList(2:end-1)), 'Color', [1 1 1]*0.4, 'LineWidth', 0.8)

    % figure setting
    set( gca, 'LineWidth', 1, 'FontSize', 14)
    xlabel( 'log(Intensity) (a.u.)'), ylabel( 'Counts')
    title( 'Signal Intensity')
    % xlim( [3 10]) % SK830


for k = 1: numel( binList)- 1

    c = c + 1;
    % ~~~~~~~~ Condition ~~~~~~~~~
    condAmp = binData >= binList(k) & binData < binList(k+1); % amp in this bin

    % xNorm = abs( xxNorm( cond));
    % lNorm = llNorm( cond);
    xNorm = abs( xNorm40( condAmp,:));
    lNorm = lNorm40( condAmp,:);

    % legtxt = sprintf( '%s: %g%%-%g%%', binName, binPer(k)*100, binPer(k+1)*100);
    legtxt = sprintf( '%g%%', binPer(k+1)*100);
    fprintf( '      %d/%d tracks,   %s\n', sum( condAmp), sum( cond), legtxt)


    % 1. plot abs( xNorm)
    figure(f1), hold on
    % [~, edges] = histcounts( xxNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
    %     tmp = movmean( edges, 2);   centers = tmp( 2:end);
    % errorbar( centers, mean( mX), std( mX), 'LineWidth', 2.5, 'DisplayName', legtxt, 'Color', colorList(c,:)), hold on
    [N, edges] = histcounts( xNorm, 'BinWidth', xNormBin, 'BinLimits', [0 1], 'Normalization', 'probability');
        tmp = movmean( edges, 2);   centers = tmp( 2:end);
    plot( centers, N, 'LineWidth', 2.5, 'Color', colorList(c,:), 'DisplayName', legtxt)

end


%% figure setting
titletxt = sprintf( '%s, %s%s', strain, strainName, extraName);

% xNorm
figure(f1)
set( gca, 'FontName', 'Arial', 'FontSize', 20)
set( gca, 'Xtick', 0:0.5:1, 'LineWidth', 1)
xlabel([ '|xNorm|, bin=' num2str( xNormBin)]), ylabel( 'Probability')
% title( titletxt)
legend( 'Location', 'northeast', 'box', 'off', 'FontSize', 11)
ylim( [0 0.15])

figure( f1)
set( gca, 'FontName', 'Arial', 'LineWidth', 1, 'FontSize', 20)
xlabel('|xNorm|'), ylabel( 'Probability')

% set( gca, 'Xtick', [0.1, 1, 10])
legend( 'Location', 'best', 'box', 'off', 'FontSize', 12)
% title( sprintf( 'EA-MSD (%d+frame, %s)', minTL, fitTxt_nl), 'FontSize', 14)



%% Save Figure to local folder

outPath = 'C:\Users\yuhuanw2\Downloads\Plots';
saveas( f1, fullfile(outPath, 'Fig2 sub_spotNorm size.svg'));