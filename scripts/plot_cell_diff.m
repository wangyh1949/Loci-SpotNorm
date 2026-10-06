%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 6/11/2026
    Last update date: 6/11/2026

    ~~~~~~ adapted from plot_spotNorm_diff.m ~~~~~~

Description: This script plots cell info based on Dα

---------------------------------------------------------------------------
%}

clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath


% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci*'));
plotNum = getPlotNum( lociList);

strainList = { 727 892 908};
nameList = { 'araC' 'Δhns' 'Δhup'};

cc = 0;
for j = plotNum
    
    cc = cc + 1;
    % set up figures for each strain
    f1 = figure( 'Position', [400+20*cc 500+20*cc 400 370]); % for xNorm
    f2 = figure( 'Position', [810+20*cc 500+20*cc 400 370]); % for LNorm
    f3 = figure( 'Position', [1220+20*cc 500+20*cc 400 370]);
    
    % load lociPos file
    load( fullfile( lociList( j).folder, lociList( j).name))

    % find strain name
    index = find( strcmp( string( strainList), strain(3:end)));    
    if isempty( index), strainName = strain;
    else, strainName = nameList{ index}; end
    

    legtxt = sprintf( '%s %s%s', strain, strainName, extraName);

    binData = Dalpha( tracksLength >= 12);    binName = 'Dα';

    % binPer = [0 0.1 0.3 0.5 0.7 0.9 1]; % binning percentage
    binPer = [0 0.05 0.1 0.4 0.8 1]; % SK830    
    binPer = [0 0.05 0.1 0.15 0.25 0.4 0.6 0.8 1]; % SK830   
    
    binList = quantile( binData, binPer); 
    colorList = flip( winter( numel( binPer))); c = 0; 
    labelName = {};

        % % display the cell numbers & spot numbers
        % fprintf( '  ~~~ %3d images, %4d cells, %5d/%5d tracks,  %5d spots,  <xNorm> = %.3f    %s-%s%s\n',...
        %     size( cellRecord, 1), sum( cellSpots > 0), sum( cond), nTracks, numel( xxNorm),...
        %     mean( xxNorm(:), 'omitnan'), expDate, strain, extraName)

    % 1. plot intensity distribution
    figure( f3)
    histogram( log( binData), 50, ... 'binWidth', 0.01,
        'EdgeColor', 'none', 'FaceAlpha', 0.5)
    xline( log( binList(2:end-1)), 'Color', [1 1 1]*0.4, 'LineWidth', 0.8)

        % figure setting
        set( gca, 'LineWidth', 1, 'FontSize', 16)
        xlabel( 'log(Dα)'), ylabel( 'Counts'), box off
        % title( 'Diffusion Coefficient')
        

    for k = 1: numel( binList)- 1

        c = c + 1;
        % ~~~~~~~~ Condition ~~~~~~~~~
        cond = binData >= binList(k) & binData < binList(k+1); % amp in this bin
        
        cInfo = cellInfo( cellNum( cond));
        nCell = numel( cInfo);

        wid = [ cInfo.width]'* 1e6; % unit: um
        length = [ cInfo.length]'* 1e6; % unit: um

        % legtxt = sprintf( '%s: %g%%-%g%%', binName, binPer(k)*100, binPer(k+1)*100);
        bintxt = sprintf( '%s: %g%%', binName, binPer(k+1)*100);
        fprintf( '      %d/%d tracks,   %s\n', sum( cond), sum( cond), bintxt)

        figure( f1)
        boxchart( c*ones( nCell, 1), wid, 'LineWidth', 1.5, 'MarkerStyle', '.'), hold on

        figure( f2)
        boxchart( c*ones( nCell, 1), length, 'LineWidth', 1.5, 'MarkerStyle', '.'), hold on
        
        labelName{ c} = sprintf( '%g%%', binPer(k+1)*100);
    end

end

%% Figure setting

% cell width
figure( f1)
set( gca, 'FontName', 'Arial', 'LineWidth', 1, 'FontSize', 16)
xlabel( sprintf( '%s Bin', binName))
ylabel( 'Cell Width (µm)'), box off
xticks( 1:c), xlim( [0.2 c+0.8])
% xticklabels( labelName)


% cell Length
figure( f2)
set( gca, 'FontName', 'Arial', 'LineWidth', 1, 'FontSize', 16)
xlabel( sprintf( '%s Bin', binName))
ylabel('Cell Length (µm)'), box off
xticks( 1:c), xlim( [0.2 c+0.8])
% xticklabels( labelName)