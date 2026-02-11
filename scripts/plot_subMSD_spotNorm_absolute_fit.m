%{
---------------------------------------------------------------------------
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 2/5/2026
    Last update date: 2/5/2026

~~~~~~ adapted from plot_subMSD_spotNorm_absolute.m & plot_subMSD_amp_fit.m ~~~~~~

Description: This script plot loci subpopulation's MSD binned by spotNorm
(xNorm or LNorm, aboslute value not percentage)

It makes 2 plots with binned data
    1. EA-MSD       (no fit)
    2. EATA-MSD     (linear fit, w/ reference line)

this script is specifically for loci tracking experiments with cell images
---------------------------------------------------------------------------
%}


clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';
lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath

% find which files to plot
lociList = dir( fullfile( lociPath, 'Loci oufti*'));
plotNum = getPlotNum( lociList);

% set strain name
strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'Ori' 'Right' 'Left' 'LacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};

minTL = 90; % minimal trackLength for plotting
fitR = 1:20;        fitTxt = sprintf( '%d:%d fit', min( fitR), max( fitR));

cc = 0;

for j = plotNum

    cc = cc + 1;

    % load loci file
    load( fullfile( lociList( j).folder, lociList( j).name))
    
    % find strain name
    index = find( strcmp( string( strainList), strain(3:end)));
    if isempty( index), strainName = strain;
    else, strainName = nameList{ index}; end
    

    % ~~~~~~~~ Condition ~~~~~~~~~
    tl = tracksLength;    frame = tracksFrame;
    cond = ( tl>= minTL) & frame(:,1) == 1;
    % cond = ( tl>= minTL);
    
    fprintf( '   plotting tracks with %d+ frames,  %d/%d tracks    %s\n\n', ...
        minTL, sum( cond), nTracks, lociName)
    
        maxTau = size( EnsMSD, 2);
        time = (1: maxTau)* timeStep;

    % divide data by spotNorm value
    binData = abs( mean( tracksxNorm40, 2));  binName = '|xNorm|'; % average of first 40 frames
    % binData = abs( mean( tracksLNorm40, 2)-0.5);  binName = 'LNorm';
    % binData = abs( tracksxNorm(:,1));  binName = '|xNorm|'; % xNorm of first spot

    binList = [0 0.1 0.2 0.3 0.4 1]; % binning absolute value
    % binList = [0 0.4 0.5 1]; % binning absolute value
    
    colorList = flip( winter( numel( binList))); c = 0;
    
    alphaFitList = nan( numel( binList)- 1, 1);
    DFitList = nan( numel( binList)- 1, 1);

    % create figures
    % f1 = figure; set( f1, 'Position', [410*cc-300 350 400 270])
    f2 = figure; set( f2, 'Position', [410*cc-300 710 400 270])

    for k = 1: numel( binList)- 1
    
        c = c + 1;

        % ~~~~~~~~ Condition ~~~~~~~~~
        condspotNorm = binData >= binList(k) & binData < binList(k+1); % spotNorm in this bin    
        condTracks = cond & condspotNorm; % overall flag for good tracks

        legtxt = sprintf( '%s: %.1f-%.1f', binName, binList(k), binList(k+1));
        % legtxt = sprintf( '%s: %.1f-%.1f, %d tracks', binName, binList(k), binList(k+1), sum( condTracks));

        fprintf( '      %d/%d tracks,   %s', sum( condTracks), sum( cond), legtxt)
        
            eaMSD = mean( EnsMSD( condTracks,:), 1, 'omitnan'); % unit: um^2
            eataMSD = mean( EnsTAMSD( condTracks, :), 1, 'omitnan');    
            

        % plotRange = 1: maxTau;
        plotRange = 1: minTL;
        
        % % 1. plot EA-MSD
        % figure( f1)
        % scatter( time( plotRange), eaMSD( plotRange), 20, 'filled', 'MarkerFaceColor', colorList(c,:), ...
        %     'DisplayName', legtxt), hold on

        % 2. plot EATA-MSD
        figure( f2)
        scatter( time( plotRange), eataMSD( plotRange), 20, 'filled', 'MarkerFaceColor', colorList(c,:), ...
            'DisplayName', legtxt), hold on

            % linear fit
            f = polyfit( log( time( fitR)), log( eataMSD( fitR)), 1);
            alphaFit = f(1);    DFit = exp( f(2))/4;    
            alphaFitList(c) = alphaFit;    DFitList(c) = DFit;
            t = time( fitR);   MSDFit =  4* DFit* t.^ alphaFit;
            fprintf( '       EATA fit: D = %.1e,  alpha = %.2f\n', DFit, alphaFit) % print fit result
    end


    % Figure Setting    
    [limX, limY] = findLim( timeStep, strain);
        % limY = 'auto';
    
    % % 1. EA-MSD
    % figure( f1), set( gca, 'LineWidth', 1, 'FontSize', 14)
    % xlabel( 'Time (s)'), ylabel( 'EA-MSD (µm^2)')
    % legend( 'Location', 'best', 'box', 'off', 'FontSize', 8)
    % title( sprintf( '%s %s (%s, %d+f)', strain, strainName, expDate, minTL), 'FontSize', 14)
    % set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
    % xlim( limX), ylim( limY)
    
    
    % 2. EATA-MSD
    figure( f2), set( gca, 'LineWidth', 1, 'FontSize', 14)
    xlabel( 'Time (s)'), ylabel( 'EATA-MSD (µm^2)')
    legend( 'Location', 'northwest', 'box', 'off', 'FontSize', 8)
    title( sprintf( '%s %s (%s)', strain, strainName, expDate), 'FontSize', 14)
    set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
    xlim( limX), ylim( limY)
    
        % set up & plot reference line
        tRef = time( fitR);
        alphaRef = alphaFitList(1);    textRef = sprintf( '\\alpha=%.2g', alphaRef);
        alphaRef2 = alphaFitList(end);    textRef2 = sprintf( '\\alpha=%.2g', alphaRef2);
        MSDRef1 =  4* DFitList(1)* 1.2* tRef.^ alphaRef;
        MSDRef2 =  4* DFitList(end)* 0.8* tRef.^ alphaRef2;
        plot( tRef, MSDRef1, 'k--', 'LineWidth', 1, 'color', [1 1 1]*0.3, 'HandleVisibility','off')
        plot( tRef, MSDRef2, 'k--', 'LineWidth', 1, 'color', [1 1 1]*0.3, 'HandleVisibility','off')
        text( tRef(end)*0.8, MSDRef1(end)*1.2, textRef, 'FontSize', 10, 'Color', [1 1 1]*0.3)
        text( tRef(end)*0.8, MSDRef2(end)*0.8, textRef2, 'FontSize', 10, 'Color', [1 1 1]*0.3)

end



%% Function

function plotNum = getPlotNum( list)

    plotNum = 1; % for multiple file of the same strain    
    
    if length( list) > 1
        fprintf( '\n~~~~~ There are %d files for this strain ~~~~~\n', length(list));
        for k = 1:length( list)
            disp( ['  ' num2str(k) '. ' list(k).name])
        end
        tmp = string( input( '\nWhich files do you want to plot (like: 1 3 4): ', 's')); % just input: 1 3 4 5
        % plotNum = double( regexp( tmp, '\d+', 'match')); % find all numbers but not extra space
        plotNum = sscanf( tmp, '%d')'; % parses any valid number format        
    end
    
fprintf( '\n')
end


function [limX, limY] = findLim( timeStep, strain)

    if timeStep == 1
        limX = [0.8 100]; limY = [5e-3 0.1]; % 100-1000ms
    elseif timeStep == 0.2
        limX = [0.1 20]; limY = [4e-3 0.06]; % 20-200ms
        if strcmp( strain, 'SK311')
            limX = [0.1 30]; limY = [1e-3 0.1]; % 20-200ms, with SK311
        end
    elseif timeStep == 0.02
        limX = [0.01 3]; limY = [3e-4 0.03]; % 20-200ms, with SK311
    else
        limX = 'auto'; limY = 'auto';
    end
end