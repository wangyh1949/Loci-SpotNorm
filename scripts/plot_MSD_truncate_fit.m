%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 1/26/2026
    Last update date: 1/26/2026

    ~~~~~~~ adapted from plot_MSDxy_compare.m & plot_MSD_fit.m ~~~~~~~

Description: This code re-calculates the TA-MSD for tracks that are
truncated to the same length (minTL).

It plots the MSD for both the original tracks and the truncated tracks in
the same figure.

I used bootstrap to estimate the uncertainty of the fitting parameters alpha
---------------------------------------------------------------------------
%}

%% 1. 

clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';

lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
tfPath   = fullfile( lociPath, 'tracksFinal'); % subfolder under lociPath

% find which files to plot
tfFile = dir( fullfile( tfPath, 'tf oufti*'));
plotNum = getPlotNum( tfFile);


% set up strain name info
strainList = { 727 728 729 730 731 734 701 662 311 725};
nameList = { 'araC' 'Ter' 'ori' 'right' 'left' 'lacZ'... %'12tetO@LacZ' ...
    '6tetO@LacY' '6tetO@lacZ' '140tetO pJZ133' '140tetO@lacZ'};


minTL = 40; % minimal trackLength for MSD calculation
fitR = 1:20;    fitTxt = sprintf( '%d:%d fit', min( fitR), max( fitR));
c = 0;

for i = plotNum

    c = c + 1;
    % load the trackFinal file
    load( fullfile( tfPath, tfFile(i).name))

    % load the corresponding loci file
    lociName = [ 'Loci oufti ' strain ' ' expDate extraName];
    load( fullfile( lociPath, lociName))
    
        % find strain name
        index = find( strcmp( string( strainList), strain(3:end)));    
        if isempty( index), strainName = strain;
        else, strainName = nameList{ index}; end

        
    nTracks = numel( tracksFinal);

    % ~~~~~~~~ Condition ~~~~~~~~~
    tl = cellfun( @length, {tracksFinal.amp}');
    cond = ( tl >= minTL) & tracksFrame(:,1) == 1;
    % cond = ( tl >= minTL);

        fprintf( '   plotting %d+ frames,  %d/%d tracks     [%s]\n', ...
            minTL, sum( cond), nTracks, tfName)

    % exclude tracks that doesn't meet the condition
    tf = tracksFinal( cond);

    % initialization for MSD_trunc
    maxT = minTL;
    EnsMSD_trunc = nan( numel( tf), maxT-1);     EnsTAMSD_trunc = nan( numel( tf), maxT-1);

    for j = 1: numel( tf)

        nFrames = numel( tf(j).amp); % number of frames in this track
        traj = tf(j).traj* pixelSize; % unit: m

        % calculate MSD for tracks with truncation (maxT)
        [ EnsMSD_trunc(j,:), EnsTAMSD_trunc(j,:)] = getMSD( traj(1:maxT,:), maxT);
    end
    
    % convert unit from m --> um
    EnsMSD_trunc = EnsMSD_trunc* 1e12;      % unit: um^2
    EnsTAMSD_trunc = EnsTAMSD_trunc* 1e12;  % unit: um^2


    % calculate EA-MSD & EATA-MSD
    eaMSD = mean( EnsMSD( cond,:), 1, 'omitnan'); % unit: um^2
    eataMSD = mean( EnsTAMSD( cond, :), 1, 'omitnan');

    eaMSD_trunc = mean( EnsMSD_trunc, 1, 'omitnan'); % unit: um^2
    eataMSD_trunc = mean( EnsTAMSD_trunc, 1, 'omitnan');


    % determine frame time (timeStep) for the movie
    frameT = str2double( erase( regexp( extraName, '\d+ms', 'match', 'once'), 'ms'));
    timeStep = frameT* 1e-3; % unit: s
    time = (1: maxT-1)* timeStep;
    plotRange = 1: maxT-1;


    % set legend text
    legtxt = sprintf( '%s, %s', strain, expDate); 
    legtxt2 = sprintf( '%s,%s', strainName, extraName);

    % set up figures
    figure, set( gcf, "Position", [410*c-110 550 400 270])
    colorList = get( gca,'colororder');
    

    % 1. plot EATA-MSD & linear fit
    scatter( time( plotRange), eataMSD( plotRange), 20, 'filled', 'MarkerFaceColor', colorList(1,:),...
        'MarkerFaceAlpha', 0.6, 'HandleVisibility','off'), hold on        
    
        [Dapp, alphaFit] = plotFitLinear( time, eataMSD, fitR, colorList(1,:), 'MSD');
        fprintf( '        MSD: Dapp= %.2e, alpha= %.2f\n', Dapp,  alphaFit);
        
        % Bootstrap
        EnsMSD_good = EnsTAMSD( cond, :); nboot = 1000;
        fit_boot = bootstrapFit( EnsMSD_good, nboot, time, fitR);
        D_boot = fit_boot(:,1);     alpha_boot = fit_boot(:,2);

    % 2. plot EATA-MSD (truncated) & linear fit
    scatter( time( plotRange), eataMSD_trunc( plotRange), 20, 'filled', 'MarkerFaceColor', colorList(2,:),...
        'MarkerFaceAlpha', 0.6, 'HandleVisibility','off')

        [Dapp_trunc, alphaFit_trunc] = plotFitLinear( time, eataMSD_trunc, fitR, colorList(2,:), 'trunc');
        fprintf( '        trunc: Dapp= %.2e, alpha= %.2f\n', Dapp_trunc,  alphaFit_trunc);

        % Bootstrap
        EnsMSD_trunc_good = EnsTAMSD_trunc; nboot = 1000;
        fit_boot_trunc = bootstrapFit( EnsMSD_trunc_good, nboot, time, fitR);
        D_boot_trunc = fit_boot_trunc(:,1);     alpha_boot_trunc = fit_boot_trunc(:,2);

    
    % figure setting
    figure( gcf), set( gca, 'LineWidth', 1, 'FontSize', 14)
    xlabel( 'Time (s)'), ylabel( 'EATA-MSD (µm^2)')
    legend( 'Location', 'southeast', 'box', 'off', 'FontSize', 12)
    title( sprintf( '%s %s (%s, %d+f)', strain, strainName, expDate, minTL), 'FontSize', 14)
    set( gca, 'Xscale', 'log', 'YScale', 'log'), box on
    [limX, limY] = findLim( timeStep, strain);
    xlim( limX), ylim( limY)

    
    % 3. plot bootstrap results for alpha
    figure, set( gcf, "Position", [410*c 190 200 270])
    boxchart( ones( nboot, 1), alpha_boot, 'MarkerStyle', '.'), hold on
    boxchart( 2*ones( nboot, 1), alpha_boot_trunc, 'MarkerStyle', '.')
    set( gca, 'XTick', [1 2], 'XTickLabel', {'reg', 'trc'}, 'FontSize', 14)
    ylabel( '\alpha from Bootstrap'), grid on
    title( sprintf( '%s %s', strain, strainName), 'FontSize', 14)
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
        plotNum = sscanf( tmp, '%d')'; % parses any valid number format
    end
    
fprintf( '\n')
end


function [ EnsMSD, EnsTAMSD] = getMSD( traj, maxT)

    nFrames = size( traj, 1);
    
    EnsMSD = nan(1, maxT-1);
    EnsTAMSD = nan(1, maxT-1);

    % ~~~~ Ensemble-Averaged MSD ~~~~~
    dr = traj( 2: min( nFrames, maxT), :) - traj( 1, :);
    EnsMSD( 1: size( dr, 1)) = sum( dr.^2, 2); % ensemble MSD for EA-MSD plotting        

    % ~~~~ Time-Averaged MSD ~~~~~
    for tau = 1: min( nFrames, maxT)-1
        dr = traj( tau+1:end, :) - traj( 1:end-tau, :);
        dr2 = sum( dr.^2, 2);
        EnsTAMSD(tau) = mean( dr2, 'omitnan'); % ensemble TA-MSD for EATA-MSD plotting
    end    
end


function [Dapp, alphaFit] = plotFitLinear( time, eataMSD, fitR, lineColor, axisName)

        % linear fit
        f = polyfit( log( time( fitR)), log( eataMSD( fitR)), 1);
        alphaFit = f(1);    Dapp = exp( f(2))/4; 
        t = time( fitR);   MSDFit =  4* Dapp* t.^ alphaFit;

        % plot fitting line
        plot( t, MSDFit, 'LineWidth', 1, 'color', lineColor, ...
            'DisplayName', sprintf( '%s: D\\alpha=%.1e, \\alpha=%.2f', axisName, Dapp, alphaFit))
end


function [limX, limY] = findLim( timeStep, strain)

    if timeStep == 1
        limX = [0.8 100]; limY = [5e-3 0.1]; % 100-1000ms
    elseif timeStep == 0.2
        limX = [0.1 20]; limY = [3e-3 0.04]; % 20-200ms
        if strcmp( strain, 'SK311')
            limX = [0.1 30]; limY = [1e-3 0.1]; % 20-200ms, with SK311
        end
    elseif timeStep == 0.02
        limX = [0.01 3]; limY = [3e-4 0.03]; % 20-200ms, with SK311
    else
        limX = 'auto'; limY = 'auto';
    end
end


function fit_boot = bootstrapFit( EnsMSD, nboot, time, fitR)

    fit_boot = nan( nboot, 2);
    n = size( EnsMSD, 1);

    tic
    for b = 1: nboot
        % Sample trajectories WITH replacement
        idx = randi( n, n, 1);
        MSD_boot = mean( EnsMSD(idx, fitR), 1, 'omitnan');

        % linear fit
        f = polyfit( log( time( fitR)), log( MSD_boot( fitR)), 1);
        alphaFit = f(1);    DFit = exp( f(2))/4;

        fit_boot(b,:) = [ DFit, alphaFit];  % store fitting result
    end
    % toc
end