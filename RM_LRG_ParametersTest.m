%% Test localization parameters LRG Super-Res
%RMN 4/29/2024

%Use LRG_SuperRes to identify particles
%test changing user defined variables
%preview localization results

clear
close all
clc

numFiles = 2; %how many files to test

fileNames = "";
filePaths = "";
for i=1:numFiles
    [file, path] = uigetfile('.tif','Select a TIF image file')
    fileNames = [fileNames, file];
    filePaths = [filePaths, path];
end



%% User Defined Parameters - Run section for 2D and 3D
%Manually tweak these parameters to test ideal fitting for your data
% General Parameters
e.runsuperres='true'; % (true or false) false = find kinetics using diffraction limited data. true = find kinetics using superlocalized data
e.startframe=100; % The first frame to analyze 
e.stopframe= 150;% The last frame to analyze
e.nframes=e.stopframe-e.startframe+1;
e.xmin=1;
e.xmax=512;
e.ymin=1;
e.ymax=512;
e.pixelSize=160; % nm

% Data Parameters
e.filename='yaxis'; % Name of the file to read in if reading in a single file. Do not include extension. 
e.loadsif='true'; % true = load data from a .sif file (if this is the first analysis for instance) false = data will be loaded from a .mat file

% Particle Identification Parameters
e.BackgroundThreshold = 1; %multiplicative minimum treshold over local background for particle ID (before adding sigma)
e.selectROI = 'false'; %prompts uuser to draw an ROI for analyzis for each file
e.wasCut = 'false'; %ROI cut was used before

e.SNR_enhance=true; % Boost the signal to noise ratio. ({true} or false)
e.local_thd= true; % Use local background. ({true} or false)
e.sigma=2; % how many stdev to add up as a threshold (Default = 3)-  higher is less lenient
e.wide2=2; % pixels Local maximum cutoff distance (Real # between 1 and 5. Default = 3) - higher is less lenient
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
e.nevent=5;% minimum peak intensity of the spot - to distinguish specific and non-specific intrs (5)
e.SRCorrFactor=0.6; %the minimum cross correlation factor between a spot and the standard Gaussian peak (0.6)
e.FinalLocatThresh=3;%how many FinalLocatSigma are used to rule out spfc vs non-spfc (3)

% Parameters for kinetics analysis
e.kinetiC = 'false'; %whether or not to calculate kinetics
e.dataSpace=0; %ms Dead time of the detector
e.dataTime=30; %ms Integration time
e.datafreq=e.dataSpace+e.dataTime; %ms frame rate
e.chngpt=0; %calculate kinetics by change point (=1) or counting molec. ID (=0)?
e.CPstate=0; % program ID change point states (0) or # for user defined states (>0)

%% Run LRG_SuperRes and Plot Image Data + Localizations
close all
clc

addpath('Z:\RicardoMongeN\MATLAB\LRG_SuperRes_Kinetics_Final\LRG_SuperRes_Kinetics_Final')
addpath('Z:\RicardoMongeN\MATLAB\LRG_SuperRes_Kinetics_Final\LRG_SuperRes_Kinetics_Final\RMN Codes')

figure(1)
for i=2:numel(fileNames); 
    e.filename=fileNames(i);
    
    e.codedir= filePaths(i);
    e.path= filePaths(i);
    addpath(e.path); 
    
    %get file size info, then limit to whole file if needed
    fileInfo = imfinfo(strcat(char(e.path),char(e.filename)));
    fileInfo = [fileInfo(1).Width(1),fileInfo(1).Height(1),size(fileInfo,1)];
%     e.stopframe=fileInfo(3);
%     e.nframes=e.stopframe-e.startframe+1;
    e.xmax=fileInfo(1);
    e.ymax=fileInfo(2);
%run localization code with chosen variables
    RM_LRG_SuperRes_Run(e);  %comment out if not running localization code,
    %and simply load in the files
    
    load(strcat(e.filename,'_analyzed.mat'));
    
    figure(1)
    y=figure(1);
        
    for k=1:ceil(size(LocatStore,2));
        figure(1)
        hold on
       imagesc(Data1(:,:,k)) %used for per frame
%         imcontrast
        colorbar%('West','Color','w')
        colormap hot
        axis image      
        if k==1 %define desired image contrast
            imcontrast
            display('Set contrast and Click on image to continue')
            waitforbuttonpress 
        end
        if k>1 %apply contrast to all subsequent frames
            caxis('manual')
        end
  
        %localizations plot - adds location markers over raw data       
        plot(LocatStore(1,k).PSFfinal(:,2),LocatStore(1,k).PSFfinal(:,1),'bo',...
            'MarkerSize',5,'LineWidth',2); %plot all particle localizations per frame
                 
        frame = getframe(y);
        im = frame2im(frame);  

%         Generate Name for locs plot(s) movie
        pltName = erase(e.filename,'.tif');
        pltName = strcat(path,pltName,'__',...
                num2str(e.sigma),'sg',...
                num2str(e.wide2),'wd',...
                num2str(e.Gauss_width),'gw',...
                num2str(e.BackgroundThreshold),'BT',...
                '_',e.fitting);
        toWrite = char(strcat(pltName,'_Movie.tif')); %will save the movie in the same folder as chosen data .tiff file
        
        if k == 1 %initial frame

          imwrite(im, toWrite)%, 'Loopcount',inf); 
        else  %attach the rest of the frames
          imwrite(im,toWrite,'WriteMode','append'); 
        end 
            
        drawnow %preview current frame -commment out for faster loading/saving
    end   
%   
% %     saveas(gcf,strcat(pltName,'_Locs.fig'));
% %    drawnow
end 

disp('Done with All Files')