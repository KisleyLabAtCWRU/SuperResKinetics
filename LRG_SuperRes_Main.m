%% LRG_SuperRes_Main

filelist=[77];

%% User Defined Parameters
% General Parameters
e.codedir='U:\Lydia\MATLAB\SuperPosition-20130215\LRG_SuperRes_Kinetics_Final\';
e.runsuperres='true'; % (true or false) false = find kinetics using diffraction limited data. true = find kinetics using superlocalized data
e.startframe=1; % The first frame to analyze
e.stopframe=5000; % The last frame to analyze
e.nframes=e.stopframe-e.startframe+1;
e.xmin=1;
e.xmax=512;
e.ymin=1;
e.ymax=512;
e.pixelSize=64; % nm

% Data Parameters
e.path='U:\Lydia\Data\SOFI121014\';
e.filename='4'; % Name of the file to read in if reading in a single file. Do not include extension. 
e.loadsif='false'; % true = load data from a .sif file (if this is the first analysis for instance) false = data will be loaded from a .mat file

% Particle Identification Parameters
e.SNR_enhance=true; % Boost the signal to noise ratio. ({true} or false)
e.local_thd=true; % Use local background. ({true} or false)
e.sigma=3; % how many std to add up as a threshold (Default = 3)
e.wide2=3; % pixels Local maximum cutoff distance (Real # between 1 and 5. Default = 3)
e.Gauss_width=2; % pixels Width of the PSF Gaussian. (real # beween 1 and 3. Default is 2)
e.fitting='rc'; % Fitting method to use. ({'rc'} 'gs' 'el') 'rc' = radial symmetry, 'gs' = gaussian fitting, 'el' = Euler fitting
e.test=false; % Generate a figure showing identified partilce locations in the first frame ('true' or {'false'}) Careful, this will be done for each frame. 

% Diffraction Limited Event Grouping Parameters
e.sameLocation=5; % pixels Distance threshhold for grouping particles in diffraction limited data (Real # Default is 2)

% SuperResolution Parameters
e.nzoom=6; % Super resolution zoom value
e.sigmarad=25; %nm sigma radius used to generate superresolution data. 

% SuperResolution Event Grouping Parameters
e.FinalLocatSigma=0.5;
e.nevent=2;
e.SRCorrFactor=0.6;
e.FinalLocatThresh=6.25;

% Parameters for kinetics analysis
e.dataSpace=0; %ms Dead time of the detector
e.dataTime=30; %ms Integration time
e.datafreq=e.dataSpace+e.dataTime; %ms frame rate
e.chngpt=0; %calculate kinetics by change point (=1) or counting molec. ID (=0)?
e.CPstate=0; % program ID change point states (0) or # for user defined states (>0)

%% Run the Analysis
for i=1:numel(filelist)
    e.filename=num2str(filelist(i));
    LRG_SuperRes_Run(e)
end
