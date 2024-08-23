%%  RM_ LRG Super Res Main Run

%Use LRG_SuperRes to identify particles over many
%frames (bining)
clear
close all
clc


startloc='Z:\RicardoMongeN\IX-73 collections\2024_03_21\s4_wKo1_0um_A1_100nM_DOX_5uL_g1_532_2.7v_30ms_400g_1\';%folder where file explorer opens
numFiles=1; %how many files you want to open
multiSelect = 'false'; %set to true if selecting multiple files in the exact same folder (no subfolders)

cd(startloc);

fileNames = "";% stuctures to store files information
filePaths = "";

if strcmp(multiSelect,'true')==0
    for i=1:numFiles
        [file, path] = uigetfile('.tif')
        
        fileNames = [fileNames, file];
        filePaths = [filePaths, path];
    end
else
   [file,path] = uigetfile(startloc,...
   'Select One or More Files', ...
   'MultiSelect', 'on');
    numFiles = size(file,2);
    for i=1:numFiles
        fileNames = [fileNames,file(i)];
        filePaths = [filePaths, path];
    end
end


disp('All Files Selected')
  %%   User Defined Parameters - Run section for 2D and 3D
% General Parameters
%e.codedir='C:\Users\Kisleylab\Desktop\UsersTemp\RicardoMN\2022_06_02\'; 
e.runsuperres='true'; % (true or false) false = find kinetics using diffraction limited data. true = find kinetics using superlocalized data

%Frames to analyze
allFrames = 'true';%set to 'false' if want to manually define range below
e.startframe=1; % The first frame to analyze 
e.stopframe=2000;% The last frame to analyze
e.nframes=e.stopframe-e.startframe+1;

%ROI size variables
fullArea = 'true'; %set to 'false' if want to manually define area below
e.xmin=1;
e.xmax=512;
e.ymin=1;
e.ymax=512;


e.pixelSize=160; % nm, corresponds to camera pixel X magnification

% Data Parameters
%e.path='Z:\RicardoMongeN\IX-73 collections\2023_03_17\';
e.filename='yaxis'; % Name of the file to read in if reading in a single file. Do not include extension. 
e.loadsif='true'; % true = load data from a .sif file (if this is the first analysis for instance) false = data will be loaded from a .mat file

% Particle Identification Parameters
e.BackgroundThreshold = 1; %multiplicative minimum treshold over local background for particle ID (before adding sigma)
e.selectROI = 'false'; %prompts uuser to draw an ROI for each file
e.wasCut = 'false'; %if data was previously cut down by ROI select

e.SNR_enhance=true; % Boost the signal to noise ratio. ({true} or false)
e.local_thd= true; % Use local background. ({true} or false)
e.sigma=2; % how many std to add up as a threshold (Default = 3)-  higher is less lenient
e.wide2=1.8; % pixels Local maximum cutoff distance (Real # between 1 and 5. Default = 3) - higher is less lenient
e.Gauss_width=2; % pixels Width of the PSF Gaussian. (real # beween 1 and 3. Default is 2) - higher is less lenient
e.fitting='rc'; % Fitting method to use. ({'rc'} 'gs' 'el') 'rc' = radial symmetry, 'gs' = gaussian fitting, 'el' = Euler fitting
e.test=false; % Generate a figure showing identified partilce locations in the first frame ('true' or {'false'}) Careful, this will be done for each frame. 

% Diffraction Limited Event Grouping Parameters
e.sameLocation=2; % pixels Distance threshhold for grouping particles in diffraction limited data (Real # Default is 2)

% SuperResolution Parameters
e.nzoom=20; % Super resolution zoom value
e.sigmarad=25; %nm sigma radius used to generate superresolution data. 

% SuperResolution Event Grouping Parameters
e.FinalLocatSigma=0.5;%the standard deviation of the single binding site - used for model 2D peak (0.5)
e.nevent=2;% minimum peak intensity of the spot - to distinguish specific and non-specific intrs (5)
e.SRCorrFactor=0.6; %the minimum cross correlation factor between a spot and the standard Gaussian peak (0.6)
e.FinalLocatThresh=3;%how many FinalLocatSigma are used to rule out spfc vs non-spfc (3)

% Parameters for kinetics analysis
e.kinetiC = 'true'; %whether or not to calculate kinetics
e.dataSpace=0; %ms Dead time of the detector
e.dataTime=30; %ms Integration time
e.datafreq=e.dataSpace+e.dataTime; %ms frame rate
e.chngpt=0; %calculate kinetics by change point (=1) or counting molec. ID (=0)?
e.CPstate=0; % program ID change point states (0) or # for user defined states (>0)

  %% Run LRG_SuperRes 
close all
clc

%change to folder locations of LRG package
addpath('Z:\RicardoMongeN\MATLAB\LRG_SuperRes_Kinetics_Final\LRG_SuperRes_Kinetics_Final')
addpath('Z:\RicardoMongeN\MATLAB\LRG_SuperRes_Kinetics_Final\LRG_SuperRes_Kinetics_Final\RMN Codes')


for i=2:numel(fileNames); %starts at 2 since first space in cell is blank
    e.filename=fileNames(i);%loads selected file information
    
    e.codedir= filePaths(i);
    e.path= filePaths(i);
    addpath(e.path); 
    
    %get file size info, then limit to whole file if needed
    %comment these lines out if want to manually assign these values in the
    %user defined parameters
    fileInfo = imfinfo(strcat(char(e.path),char(e.filename)));
    fileInfo = [fileInfo(1).Width(1),fileInfo(1).Height(1),size(fileInfo,1)];
    
    %check if user defined ranges of defaults
    if strcmp(allFrames,'true')
        e.stopframe=fileInfo(3);
        e.nframes=e.stopframe-e.startframe+1; %defaults to analyzing the whole image
    end
    if strcmp(fullArea,'true')
        e.xmax=fileInfo(1); %assumes full ROI initially
        e.ymax=fileInfo(2);
    end
%run localization analyzis
    RM_LRG_SuperRes_Run(e);  %custom version of main script
end
disp('All Done')