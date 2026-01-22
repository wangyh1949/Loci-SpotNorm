%{
Author: Yu-Huan Wang (Kim Lab at UIUC) - yuhuanw2@illinois.edu
    Creation date: 3/25/2024
    Last update date: 12/1/2025

    ~~~~~~~ adapted from loci_script.m ~~~~~~~

Description: This script prepares the raw loci tracking data to be ready
for uTrack/oufti analysis. 

~~~~~~~~~~~~~~ folder & file naming rule ~~~~~~~~~~~~~~
1. Raw data folder in hard drive 'F:\2025\250922-SK734'
        (should follow the exact folder format)
        i.e. 'expDate-strain extraName'
              extraName examples: Gly, rif, IPTG, etc
2. Inside raw folder, should contain 'tracking00x' folders with exported tiff files
3. Inside raw folder, should contain a 'phase' folder containing phase tiff files
4. Should have a local folder storing the analysis result from uTrack & oufti
        i.e. 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data' or anywhere you like
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

this script will create a folder with the same folderName of the raw data
folder, and copy all tiff images, and create empty tracking00x folders that
corresponds to the raw data folder.

The empty tracing00x folders will later be used as output folders to save
uTrack analysis results

---------------------------------------------------------------------------
%}

clear, clc, close all
 
% lociPath = uigetdir( 'E:\2025', 'Please choose the hard drive folder that contains raw tracking data');
% [~, folderName] = fileparts( lociPath); % get the raw data folder name

folderName = '260115-SK892'; % folderName for the tracking experiment
lociPath = fullfile( 'F:\2026\', folderName); % folder containing expermental data (on a hard drive)

% corresponding uTrack output folder on the computure
% change to your own local folder that stores the analysis results
dataPath = [ 'C:\Users\yuhuanw2\Documents\MATLAB\Lab Data\2026\' folderName];

tmp = dir( fullfile( lociPath, '*racking*'));
% tmp = dir( fullfile( lociPath, '*\*racking*'); % for tracking folders are in some subfolder

list = tmp( [ tmp.isdir]); % folders: tracking001, tracking002, ..., tracking00x


%% 1. Copy the BF or phase images to local folder

tifFile = dir( fullfile( lociPath, '*.tif'));

if ~isempty( tifFile)
    % move all epi00x.tif into the 'phase' subfolder
    copyfile( fullfile( lociPath, '*.tif'), fullfile( dataPath, 'phase'))
    fprintf( '\n ~~~ TIFF images copied    %s ~~~\n\n', lociPath)
elseif ~isempty( dir( fullfile( lociPath, 'ph*')))
    % copy the whole folder 'phase' 
    copyfile( fullfile( lociPath, 'ph*'), dataPath)
    fprintf( '\n ~~~ phase images copied    %s ~~~\n\n', lociPath)
else
    fprintf( '\n ~~~ No TIFF images found ~~~\n\n')
end


%% 2. Create corresponding empty tracking folders for uTrack output

for i = 1: length( list)
    newPath = fullfile( dataPath, list(i).name);
    if ~exist( newPath, 'dir')
        mkdir( newPath)
        fprintf( ' ~~~ %s folder created ~~~\n', list(i).name)
    else
        fprintf( ' ~~~ %s folder already existed ~~~\n', list(i).name)
    end
end

cd( dataPath)


% Run uTrack for loci tracking    
fprintf( '\n~~~~~~    Time to run uTrack for the tracking    ~~~~~~\n')


% for loci tracking data, I'm currently using following parameter for spot detection
%   Gaussian std = 1 pixel (Andor camera)
%   Alpha (Detection) = 0.01 
%   Alpha (fitting) = 0.01
% 
%   Tracking: 0 gap, 40+ frames, 3 pix (depends on specific strain)




