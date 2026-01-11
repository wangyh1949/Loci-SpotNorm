%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 1/11/2026
    Last update date: 1/11/2026

    ~~~~~~~ adapted from combLociTF.m & lociPos_comb.m ~~~~~~~
    ~~~~~~~ adapted from lociAnalysis_spotNorm.m ~~~~~~~

Description: This code combines multiple Loci_spotNorm single-day files and
does diffusion analysis for the combined data.

Single-day files are moved to the 'single day' subfolder (relative path)
for both tf & loci files.

Two files are combined and saved
    1. tf oufti: contains raw track information tracksFinal 
    2. Loci oufti: contains quantities calculated from diffusion analysis

---------------------------------------------------------------------------
%}

%% 1. Combine tf oufti files

clear, clc, close all

% set the path for storing analysis results
varPath = 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\Track and Cell Variables\';

lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
tfPath   = fullfile( lociPath, 'tracksFinal'); % subfolder under lociPath


% find which strain to combine
strain = getCombStrain( tfPath);

% find which files to combine
tfFile = dir( fullfile( tfPath, ['tf oufti*' strain '*']));
combNum = getCombNum( tfFile);


% initialization of combining quantities
numComb = length( combNum);

tfComb = cell( numComb, 1);     cellMesh = cell( numComb, 1);
combRec = cell( numComb, 1);    cellRec = cell( numComb, 1);
num = [0 0];    c = 0;

for i = combNum
    
    c = c + 1;
    % load the tf oufti file
    tfPath   = fullfile( lociPath, 'tracksFinal'); % subfolder under lociPath
    load( fullfile( tfPath, tfFile(i).name)),  fprintf( '    ~~~ %s loaded ~~~\n', tfName)

    lociName = [ 'Loci oufti ' strain ' ' expDate extraName];
    
    % combine the tf
    tfComb{ c} = tracksFinal;
    cellMesh{ c} = cellMeshAll;
    cellRec{ c} = cellRecord;

    nTracks = size( tracksFinal, 1);
    num = num(2) + [1 nTracks]; % record the track numbers from each movies
    
    % record all info into a single cell variable
    combRec{c} = { dataPath, folderName, tfName, lociName, sigToPhoton, num}; 
end

tf = vertcat( tfComb{:});
cellMeshAll = vertcat( cellMesh{:});
cellRecord = vertcat( cellRec{:});
combRec = vertcat( combRec{:});

% reorder the cell number
tracksFinal = reNumCells_yh( tf); 

fprintf( '\n~~~~~~~ tracksFinal Combination Finished ~~~~~~~\n');


lociPath = fullfile( varPath, 'Loci SpotNorm'); % subfolder under varPath
tfPath   = fullfile( lociPath, 'tracksFinal'); % subfolder under lociPath

% get experimental info
expDate = input( '\nWhat is the combDate ( like ''240610 comb''):  ', 's');

% save tracksFinal files
tfName = ['tf oufti ' strain ' ' expDate extraName];

save( fullfile( tfPath, tfName), 'varPath', 'lociPath', 'dataPath', 'imgPath', 'tracksFinal', ...
    'cellMeshAll', 'cellRecord', 'cameraFlag', 'pixelSize', 'sigToPhoton', ...
    'combRec', 'folderName', 'expDate', 'strain', 'extraName', 'tfPath', 'tfName')

fprintf( '\n ~~~ tracksFinal:  %s saved  under  ''tfPath'' ~~~\n\n', tfName)


%% 2. Move individual 'tf oufti' files to 'Individuals' folder

% Set destination folders for single-day files
tfDstFolder = fullfile( tfPath, 'single day');
lociDstFolder = fullfile( lociPath, 'single day');

% Create subfolder if it doesn't exist
if ~exist( tfDstFolder, 'dir'), mkdir(tfDstFolder); end
if ~exist( lociDstFolder, 'dir'), mkdir(lociDstFolder); end

files = tfFile( combNum);

for k = 1: numel( files)
    % move tf files into subfolder
    movefile( fullfile( tfPath, files(k).name), tfDstFolder)
    fprintf( '   ''%s'' moved to ''tfPath\\single day'' folder \n', files(k).name)
    
    % move loci files into subfolder
    lociFileName = [ combRec{k,4} '.mat'];
    movefile( fullfile( lociPath, lociFileName), lociDstFolder)
    fprintf( '   ''%s'' moved to ''lociPath\\single day'' folder \n\n', lociFileName)    
end


%% 3. Diffusion analysis

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
    'combRec', 'folderName', 'expDate', 'strain', 'extraName', 'tfPath', 'tfName', 'lociName')


fprintf( ' ~~~ Loci file:  %s saved under  ''lociPath'' ~~~\n\n', lociName)

% cd( lociPath)


%% 4. Display analysis result summary

fprintf( ['\n ~~~~~~ %s Analysis Result ~~~~~~\n' ...
    '    %d total cells,  %d/%d valid tracks\n'...
    '    %d tracks (1st spot),  %d tracks (first 40 spots)    not in cap region \n'],...
    folderName, totalCells, sum( ~isnan( tracksxNorm(:,1))), nTracks, ...
    sum( tracksMid(:,1)), sum( min( tracksMid40, [], 2)))
    
% disp( cellRecord)



%% Function

function strain = getCombStrain( tfPath)

    % list all available strains for combination
    files = dir( fullfile( tfPath, 'tf oufti*'));
    names = { files.name};
    tokens = regexp( names, 'SK\d+', 'match');
    strainList = [ tokens{:}];
    SKlist = unique( strainList);
    
    listN = length( SKlist); 
    if listN == 0
        error('   No strain is found in the current directory \n')
    else
        fprintf( '\n~~~~~ There are %d strains ~~~~~\n', listN)
        for k = 1: listN
            fprintf( '    -  %s\n', SKlist{k})
        end
        % find which strain to combine
        tmp = input( '  Which strain do you want to combine? (e.g. ''SK662'')  ', 's');
        strainNum = strcmp( tmp, SKlist);
        if sum( strainNum) ~= 1
            error( '   ~~~~~ No strain found for input ~~~~~\n')
        else
            strain = SKlist{ strainNum};
        end
    end
end