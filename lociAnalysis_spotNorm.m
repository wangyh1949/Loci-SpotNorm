%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu 
    Creation date: 12/22/2025 
    Last update date: 12/22/2025

~~~~~~~~~~~ adapted from lociAnalysis.m & autoSpotNorm.m ~~~~~~~~~~~
  this version integrates the oufti info, deleted drift correction


Description: This script combines the output from oufti & u-track analyses
for multiple movies. It further assigns tracks to cells and calculates
normalized positions of foci. From the track information, it does diffusion
analysis (MSD, step displacement, ...)

Two files are saved
    1. tracksFinal: raw track information
    2. Loci matrix: quantities calculated from diffusion analysis


~~~~~~~~~~~~~~~~~~~~~~ Before Running the Code ~~~~~~~~~~~~~~~~~~~~~~
1. change current working directory to the folder that contains tracking00x
folders (uTrack output) and mesha.mat (oufti output)
2. change 'varPath' to your own local folder to store analysis results (all
results will be saved under the lociPath)
3. make sure the oufti output file is cleaned up and named as mesha.mat
---------------------------------------------------------------------------
%}


%% 1. Run spotNorm and combine multiple movies with oufti & uTrack results

clear, clc, close all

% set the path to store analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';


goodCellFlag = true; % true: assumes all cells are good;  false: manually select good cells
imgSaveFlag = false; % true: save overlay image of cells & tracks


lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
tfPath   = fullfile( lociPath, 'tracksFinal'); % subfolder under lociPath

if ~exist( lociPath, 'dir'), mkdir( lociPath), end
if ~exist( tfPath, 'dir'), mkdir( tfPath), end

dataPath = pwd; % current folder with raw uTrack & oufti output

% get experimental info from the folderName (see function below)
[ folderName, expDate, strain, extraName] = getExpInfo( dataPath);


% ~~~~~ 1. Find oufti & uTrack output ~~~~~

    % list all tracking folders in the directory
    trackList = dir( '*racking*');  listN = length( trackList); 
    if listN == 0
        error('   No tracking folder found in the current directory \n')
    else
        fprintf( '\n~~~~~ There are %d tracking files ~~~~~\n', listN)
        for k = 1: listN
            fprintf( '%3d. %s\n', k, trackList( k).name)
        end    
        combNum = input( '\nWhich tracks do you want to combine (0: all, or like: 3:7):   ');
        % combNum = 0;
        if combNum == 0, combNum = 1: listN; end
    end

    % list all phase contrast images
    phaseList = dir( fullfile( 'phase', 'epi*.tif'));
    imgPath = phaseList(1).folder;

    % load oufti mesh file
    tmp = 'mesha.mat'; % contains oufti results for all PC images
    if isfile( tmp)
        load( tmp, 'cellList', 'cellListN')
    else 
        error('  File "%s" not found in the current directory \n', tmp);
    end
    
    % create folder to save analysis folder
    if logical( imgSaveFlag)
        fprintf( '~~~ Img save flag is on, cell & track overlay image will be saved ~~~\n\n')
        if ~exist( 'Mesh & Tracks Image', 'dir')
            mkdir( 'Mesh & Tracks Image');
        end
    else
        fprintf( '~~~ Img save flag is off, cell & track overlay image will NOT be saved ~~~\n\n')
    end


% ~~~~~ 2. Run spotNorm and combine multiple movies ~~~~~

num = length( combNum);         tfComb = cell( num, 1);
cellMeshComb = cell( num, 1);    cellRecord = cell( num, 4);
c = 0;

for m = combNum

    c = c + 1;
    
    % run spotNorm on an individual track
    [ tfGood, cellMesh] = spotNorm_loci( trackList, phaseList, m, ...
        cellList, cellListN, goodCellFlag, imgSaveFlag);
          
    cellNum = [ tfGood.cellNumber]';    % cell number of each good track
    goodCells = unique( cellNum)';      % all goodCells are not empty and have tracks inside
    
    imgName = phaseList( m).name;
    fileName = [ folderName ' ' imgName(1:end-4)]; % (i.e. '251004-SK727 epi001')


    % store filename in new fields of the tfGood structure
    nTracks = numel( tfGood);           value = cell( nTracks, 1);    
    value(:) = { string( fileName)};    [tfGood.origin] = value{:};    

    
    % stack tfGood & cellMesh from different movies together
    tfComb{ c} = tfGood;
    cellMeshComb{ c} = cellMesh( goodCells);
    cellRecord( c,:) = {fileName goodCells folderName imgName};

    % fprintf( '        Included Good cells with tracks: %s\n\n', join( string( goodCells), ' '))
end

tf = vertcat( tfComb{:}); % concatenates tf vertically to a big structure vairable
cellMeshAll = vertcat( cellMeshComb{:});

fprintf( '~~~~~~~ tracksFinal Combination Finished ~~~~~~~\n\n');


% ~~~~~ 3. Clean up tracksFinal and save into tf files ~~~~~

tracksFinal = reNumCells_yh( tf); % reorder the cell number

% remove original uTrack output fields in tracksFinal to save space  
tracksFinal = rmfield( tracksFinal, {'tracksFeatIndxCG' 'tracksCoordAmpCG' 'seqOfEvents'});

% Calculate conversion factor 'sigToPhoton' based on camera gain
setCameraGain


    % add extraNote to the fileName (usually about imaging setting): 
    %    1000ms or 20-200ms (20ms exposure, 200ms interval)
    fprintf( '\n ~~~~~~ fileName:  %s  ~~~~~~\n', ['tf oufti ' strain ' ' expDate extraName])
    fprintf( ['\n     add extraNote to the fileName (usually about imaging setting)\n' ...
        '         1000ms or 20-200ms (20ms exposure, 200ms interval)\n\n'])
    extraNote = input( ' ~~~ extra notes: what is the imaging setting (i.e. 20-200ms):  ', 's');
    % extraNote = '20-200ms';
    
    if ~isempty( extraNote)
        extraName = [ extraName ' ' extraNote];
        % extraName = [ ' ' extraNote];
    end


% save tracksFinal files
tfName = ['tf oufti ' strain ' ' expDate extraName];

save( fullfile( tfPath, tfName), 'varPath', 'lociPath', 'dataPath', 'imgPath', 'tracksFinal', ...
    'cellMeshAll', 'cellRecord', 'cameraFlag', 'pixelSize', 'sigToPhoton', ...
    'folderName', 'expDate', 'strain', 'extraName', 'tfPath', 'tfName')

fprintf( '\n ~~~ tracksFinal:  %s saved  under  ''tfPath'' ~~~\n\n', tfName)



%% 2. Diffusion analysis

% ~~~~~ 1. Diffusion analysis: Calculate MSD ~~~~~

% determine frame time (timeStep) for the movie
%       100-1000ms  --->  frame time: 1000ms
frameT = str2double( erase( regexp( extraName, '\d+ms', 'match', 'once'), 'ms'));
% frameT = double( input( '   ~~~ What is the frame time (ms) :  '));
timeStep = frameT* 1e-3; % unit: s    
fprintf( '        ~~~~ frame time is  %d ms ~~~~\n', timeStep*1e3)

tl = cellfun( @length, {tracksFinal.amp}');
maxT = max( tl); % MSD truncation length ( min( nFrames, maxT))
fprintf( '   Longest track has %d frames, chose maxT = %d for MSD\n\n', max( tl), maxT)

% initialzation
nTracks = numel( tracksFinal);
EnsMSD = nan( nTracks, maxT-1);     EnsTAMSD = nan( nTracks, maxT-1);
alpha = nan( nTracks, 1);           Dalpha = nan( nTracks, 1);
posStd = nan( nTracks, 1);          tracksAmp = nan( nTracks, maxT);

origin = unique( [tracksFinal.origin])';    

for i = 1: nTracks
    
    nFrames = length( tracksFinal(i).amp); % number of frames in this track
    tracksAmp( i, 1:nFrames) = tracksFinal(i).amp'; % converted to photon
    posStd(i) = tracksFinal(i).std(1)* pixelSize; % std of position, unit: m
    
    % Stage Drift Correction
    frame = tracksFinal(i).frame(1): tracksFinal(i).frame(2);
    
    traj = tracksFinal(i).traj* pixelSize; % unit: m
    
    % % ~~~~ Jump Displacement ~~~~    unit: um
    % singleSteps = traj( 2:end, 1:2)- traj( 1:end-1, 1:2); % single steps
    % tracksFinal(i).steps = singleSteps'* 1e6; % unit: um
    
    % ~~~~ Ensemble-Averaged MSD ~~~~~
    dr = traj( 2: min( nFrames, maxT), :) - traj( 1, :);
    EnsMSD( i, 1: size( dr, 1)) = sum( dr.^2, 2); % ensemble MSD for EA-MSD plotting        

    % ~~~~ Time-Averaged MSD ~~~~~
    for tau = 1: min( nFrames, maxT)-1
        dr = traj( tau+1:end, :) - traj( 1:end-tau, :);
        dr2 = sum( dr.^2, 2);
        EnsTAMSD(i, tau) = mean( dr2, 'omitnan'); % ensemble TA-MSD for EATA-MSD plotting
    end    

    % TA-MSD fitting to get D & alpha,  MSD = 4Dt^alpha
    fitRange = 1: floor( (nFrames-1) / 4); fitTxt = '25% fit';
    % fitRange = 1:5; fitTxt = '1:5 fit';
        time = fitRange* timeStep;
        f = polyfit( log( time), log( EnsTAMSD( i, fitRange)), 1);
        alpha(i) = f(1);    Dalpha(i) = exp( f(2))/4;
end

% convert unit from m --> um
EnsMSD = EnsMSD* 1e12;      % unit: um^2
EnsTAMSD = EnsTAMSD* 1e12;  % unit: um^2
Dalpha = Dalpha* 1e12;      % unit: um^2/s
posStd = posStd* 1e9;       % unit: nm

% get matrix variables from the structure
tracksLength = cellfun( @length, {tracksFinal.amp}');
tracksOrigin = [ tracksFinal.origin]';
tracksFrame  = [ tracksFinal.frame]';
% steps        = [ tracksFinal.steps]'; % unit: um

cellNum = [ tracksFinal.ModCellNum]';
totalCells = max( cellNum);


% ~~~~~ 2. calculate cell geometry & spotNorm ~~~~~

% get cell information
[ cellInfo, poleBounds] = getCellInfo( tracksFinal, cellMeshAll, pixelSize);

    % initialization
    tracksLNorm   = nan( nTracks, 4);    tracksxNorm   = nan( nTracks, 2);
    tracksLNorm40 = nan( nTracks, 40);   tracksxNorm40 = nan( nTracks, 40);

    for i = 1: nTracks    
        LNorm = tracksFinal(i).spotPosNorm(:,1);
        xNorm = tracksFinal(i).spotPosNorm(:,2);    
        tracksLNorm(i,:) = [ LNorm(1), LNorm(end), min( LNorm), max( LNorm)];
        tracksxNorm(i,:) = [ xNorm(1), xNorm(end)]; % [first pt, end pt]        
        % min trackLength of loci tracks is 40 frame
        tracksLNorm40(i,:) = LNorm(1:40);
        tracksxNorm40(i,:) = xNorm(1:40);
    end

% cell subregion constraints (cell pole or middle part)
bound = poleBounds( cellNum);   % cell pole bound for each track
tracksMid = false( nTracks, 2); % based on [first spot, whole track]

LNorm = tracksLNorm; % [first spot, end spot, minLNorm, maxLNorm]

tracksMid(:, 1) = LNorm(:,1) > bound & LNorm(:,1) < 1-bound; % first spot not at pole
tracksMid(:, 2) = LNorm(:,3) > bound & LNorm(:,4) < 1-bound; % whole track not at pole

tracksMid40 = tracksLNorm40 > bound & tracksLNorm40 < 1-bound;

% save variables for plotting, tracksFinal not saved (save space)
lociName = [ 'Loci oufti ' strain ' ' expDate extraName];

save( fullfile( lociPath, lociName), 'varPath', 'lociPath', 'dataPath', 'imgPath', ...
    'cellRecord', 'cameraFlag', 'pixelSize', 'sigToPhoton', ...
    'maxT', 'timeStep', 'cellNum', 'totalCells', ...
    'nTracks', 'EnsMSD', 'EnsTAMSD', 'fitTxt', 'alpha', 'Dalpha', 'posStd', ...
    'tracksAmp', 'tracksOrigin', 'tracksLength', 'tracksFrame', ...% 'steps', ...    
    'cellInfo', 'poleBounds', 'tracksLNorm', 'tracksxNorm', 'tracksMid', ...
    'tracksLNorm40', 'tracksxNorm40', 'tracksMid40', ...
    'folderName', 'expDate', 'strain', 'extraName', 'tfPath', 'tfName', 'lociName')


fprintf( ' ~~~ Loci file:  %s saved under  ''lociPath'' ~~~\n\n', lociName)

% cd( lociPath)


%% 3. Display analysis result summary

fprintf( ['\n ~~~~~~ %s Analysis Result ~~~~~~\n' ...
    '    %d total cells,  %d/%d valid tracks\n'...
    '    %d tracks (1st spot),  %d tracks (first 40 spots)    not in cap region \n'],...
    folderName, totalCells, sum( ~isnan( tracksxNorm(:,1))), nTracks, ...
    sum( tracksMid(:,1)), sum( min( tracksMid40, [], 2)))
    
% disp( cellRecord)



%% Functions

function [ folderName, expDate, strain, extraName] = getExpInfo( trackPath)
% Extract experimental info from meshPath & meshName
% Naming matters! please follow the format rule example for the folder:
% '230712-SK187' or '230712-SK187 Suc'  

    tmp = split( trackPath, '\');

    twoLociFlag = any( contains( tmp, 'ch'));         
    folderName = tmp{ end - double( twoLociFlag)};

    tmp2 = split( folderName, ["-" " "]);   % split '210920-SK71 Suc'
    expDate = tmp2{1};                      % date of the experiment '210920'
    strain = tmp2{2};                       % strain number 'SK71'
    extraName = '';

    if length(tmp2)>=3
        extraName = erase( folderName, [ expDate '-' strain]); % everything afterwards
        % extraName = [ ' ' tmp2{3}];         % ' Loci' or ' Glu'
    end

    if contains( tmp{end}, 'GFP')           % two loci: 2 channels
        extraName = [ ' GFP' extraName];
    elseif contains( tmp{end}, 'mCh')
        extraName = [ ' mCh' extraName];
    end

end


% function plotDrift( drift)
% % this function plot the stage drift over time
% 
%     figure( 1)
%     set( gcf, 'Position', [ 800 400 450 450])
% 
%     % Plot drift with color code
%     x = drift(:,1)'; y = drift(:,2)'; z = zeros( size(x));
%     lineColor = 1: length(x);  % color code
% 
%     surface( [x;x], [y;y], [z;z], [lineColor; lineColor],...
%         'EdgeColor', 'interp', 'LineWidth', 2) % 'EdgeAlpha', 0.2
%     hold on
% 
%     % plot( movmean( x, [5 0]), movmean( y, [5 0]), 'r-', 'LineWidth', 1.5)        
%     set( gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
%     grid on, box on, axis equal
%     % xlim( [-0.5 3.5]), ylim( [-0.5 3.5])
%     title( 'Stage Drift', 'FontSize', 14)
% end

