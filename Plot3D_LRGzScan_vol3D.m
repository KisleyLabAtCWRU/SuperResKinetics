%% Sphere Test_V2
%Plot localizations for a long z-scan
%First, obtain Log of Z positions
%RicardoMN - KL
clear;close all;clc;

[file,path] = uigetfile('*.txt','Select Meta Data File');
addpath(path);
 
Digits = 2; %significant digits in measurements, currently to 10nm scale

fileID=fopen(path,'rt'); %open the file
A = fileread(file);

startPat ='"ZPositionUm": '; %mark Start pattern to search for in File array
endPat = ','; %mark End pattern to search for in File array

newA = extractBetween(A,startPat,endPat); %generate cell array, containing strings located between the Start and End patterns.
B = newA'; %take transpose of cell array

C = [Inf 1]; %initialize new array conta ining numerical form of B
for i =1:length(B)
    D = cell2mat(B(i)); %convert the cell arrays to ordinary arrays
    format longg
    C(i) = str2double(D);  %convert str arrays to numerical values
end
C2 = C - min(C);
R= round(C2,Digits);%'significant'); %round position array to the correct number of significant digits
disp('R is position vector, units of um')
% R = R.*1000; %convert to nm

clearvars -except R
%% Plot Localizations 
%RMN 3/11/2022

%Use LRG_SuperRes to identify particles over many
%frames (bining) - Or use and skip to 3D plotting (bottom)
%***Load ZLog file (made with "MakeZLog" script) AND 3D Localization Outputs for 3D Plot Comparisons
close all

numFiles = 1;

fileNames = "";
filePaths = "";
for i=1:numFiles
    [file, path] = uigetfile('.tif','Select a TIF image file')
    fileNames = [fileNames, file];
    filePaths = [filePaths, path];
end



%% User Defined Parameters - Run section for 2D and 3D
% General Parameters
% e.codedir='Z:\RicardoMongeN\IX-73 collections\2022_03_10\c1_0to6_0.25um_3.3V_30ms_300g_500ps'; 
e.runsuperres='true'; % (true or false) false = find kinetics using diffraction limited data. true = find kinetics using superlocalized data
% e.startframe=2000; % The first frame to analyze 
% e.stopframe=2200;% The last frame to analyze
% e.nframes=e.stopframe-e.startframe+1;
e.xmin=1;
e.xmax=512;
e.ymin=1;
e.ymax=512;
e.pixelSize=160; % nm

% Data Parameters
% e.path='C:\Users\Kisleylab\Desktop\UsersTemp\RicardoMN\2021_06_23';
% e.filename='yaxis'; % Name of the file to read in if reading in a single file. Do not include extension. 
e.loadsif='false'; % true = load data from a .sif file (if this is the first analysis for instance) false = data will be loaded from a .mat file

% Particle Identification Parameters
e.BackgroundThreshold = 1; %multiplicative minimum treshold over local background for particle ID (before adding sigma)
e.selectROI = 'false'; %prompts uuser to draw an ROI for each file
e.wasCut = 'false';

e.SNR_enhance=true; % Boost the signal to noise ratio. ({true} or false)
e.local_thd= true; % Use local background. ({true} or false)
e.sigma=3; % how many std to add up as a threshold (Default = 3)-  higher is less lenient
e.wide2=3; % pixels Local maximum cutoff distance (Real # between 1 and 5. Default = 3) - higher is less lenient
e.Gauss_width=2; % pixels Width of the PSF Gaussian. (real # beween 1 and 3. Default is 2) - higher is less lenient
e.fitting='rc'; % Fitting method to use. ({'rc'} 'gs' 'el') 'rc' = radial symmetry, 'gs' = gaussian fitting, 'el' = Euler fitting
e.test=false; % Generate a figure showing identified partilce locations in the first frame ('true' or {'false'}) Careful, this will be done for each frame. 

% Diffraction Limited Event Grouping Parameters
e.sameLocation=2; % pixels Distance threshhold for grouping particles in diffraction limited data (Real # Default is 2)

% SuperResolution Parameters
e.nzoom=20; % Super resolution zoom value
e.sigmarad=25; %nm sigma radius used to generate superresolution data. 

% SuperResolution Event Grouping Parameters
e.FinalLocatSigma=0.5;%loc accuracy-the standard deviation of the single binding site - used for model 2D peak (0.5)
e.nevent=5;%events-minimum peak intensity of the spot - to distinguish specific and non-specific intrs (5 events)
e.SRCorrFactor=0.6; %shape-the minimum cross correlation factor between a spot and the standard Gaussian peak (0.6)
e.FinalLocatThresh=3;%search radius - how many FinalLocatSigma are used to rule out spfc vs non-spfc (3 sigma)

% Parameters for kinetics analysis
e.kinetiC = 'true'; %whether or not to calculate kinetics
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
G = [[],[],[],[]];

for i=2:numel(fileNames); 
    e.filename=fileNames(i);
    
    e.codedir= filePaths(i);
    e.path= filePaths(i);
    addpath(e.path); 
    
    %e.sigma=cSG(i);%1.9; % how many std to add up as a threshold (Default = 3)-  higher is less lenient
    %e.wide2=cSG(i); % pixels Local maximum cutoff distance (Real # between 1 and 5. Default = 3) -lower is less lnt
   % e.Gauss_width=cSG(i); % pixels Width of the PSF Gaussian. (real # beween 1 and 3. Default is 2) -lower is more lnt

    %get file size info, then limit to whole file if needed
    fileInfo = imfinfo(strcat(path,file));
    fileInfo = [fileInfo(1).Width(1),fileInfo(1).Height(1),size(fileInfo,1)];
    %e.stopframe=fileInfo(3);
    e.xmax=fileInfo(1);
    e.ymax=fileInfo(2);
    
    disp('Starting')
    
    disp('Generating and Cutting Data as .mat as needed')  
    
    %find unique indeces of R(Z position vector)
    [R2,ia,ic] = unique(R(1,1:fileInfo(3))); %[unique positions, R2=R(ia),R=R2(ic)]
    
    %Process and plot data for specific ranges
    
    for j=1:(size(R2,2)-1); %j is unique Z position in list R2
        
        
        
        
        if j==1 
            %full movie range
            
            e.stopframe=fileInfo(3);
            e.startframe = 1;
            e.nframes=e.stopframe-e.startframe+1;
            
            %RM_LRG_SuperRes_Run(e);  %comment out if not running localization code,
            %and simply load in the files
            disp('Loading All Analyzed File')
%            load(strcat(e.filename,'_analyzed.mat'));
           % load(strcat(e.filename,"_",num2str(R2(round(size(R2,2)/2))),'um_analyzed.mat'));
           %Data1=TiffReadRM(e.filename,e.path,e.startframe,e.stopframe); 
           if strcmp(e.loadsif,'true')==1;
                disp('Converting .tif File')
                LRG_SuperRes_Read_Data(e) %moded to ROI cutting for loops
                disp('Converted')
                load(strcat(e.path,e.filename,'.mat'),'Data1');
           else
               load(strcat(e.path,e.filename,'.mat'),'Data1');
           end

           e.loadsif = 'false';
        end
        
        %Dynamic frame ranges
        e.startframe = ia(j);
        e.stopframe = ia(j+1)-1;%frame right before where position changed
        e.nframes=e.stopframe-e.startframe+1;
        e.kinetiC = 'true';
        e.loadsif = 'false';
        if isfile(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'))
            wasDone = 1;
        else 
            wasDone = 0;
        end
        
        if wasDone == 0;
            e.loadsif = 'false';
            disp('Analyzing current slice')
            RM_LRG_SuperRes_Run(e);  %comment out if not running localization code,
            %and simply load in the files
            disp('Loading Current Analyzed slice')
            load(strcat(e.filename,'_analyzed.mat'));
            save(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'),'LocatStore','e','Data1','GroupLocat','Ensemble'); 
            load(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'));
            
        else
            disp(strcat("Loading ",num2str(R2(j)),'um Analyzed File'))
            
            load(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'));
            e.loadsif = 'false';
        end

        Fm = e.startframe:1:e.stopframe; %frame numbers to match limits

        G(:,:,1,j)= uint16(LRG_SuperRes_GenerateSR(LocatStore,e));
        
    end
    
    
    %writecell(DwellSt,strcat((e.path),'DwellSt.csv'));
    disp('Done')
    
end

%writecell(DwellSt,strcat((e.path),'DwellSt.csv'));
save(strcat(e.path,"3DmapData.mat"),'G','R2','e','R');

clearvars -except G R R2 e
disp('Saved 3D Map Data');
%%  Vold3d plotting
close all
bWaist = 2000; %in nm
G = squeeze(G());% original dataset
R3 = R2(1:size(G,3)).*1000; % positions in nm
cvFt = (e.pixelSize/e.nzoom); %conversion factor
%G = uint16(G);
G2 = [[],[],[]];
R4 = [];
gCount = 1;
for i=1:size(R3,2)-1
    G2(:,:,gCount) = G(:,:,i);
    R4 = [R4, R3(i)];
    dR = R3(i+1)-R3(i);
    if dR> bWaist %only if beam waist is smaller than Z step size
        R4 = [R4, R3(i)+bWaist];
        gCount = gCount +1;
        G2(:,:,gCount) = zeros(size(G,1),size(G,2));%make slice of zeros from (z+bWaist) to z+deltaR
    end

    if i==size(R3,2)-1
        G2(:,:,gCount+1) = G(:,:,i+1);
        R4 = [R4, R3(i+1)];
    end
    gCount = gCount +1;
end

%%
addpath('Z:\RicardoMongeN\MATLAB')
figure()
%clim([0 25])
vol3d('Cdata',G2(:,:,20:25),...
    'ZData',R4(1,20:25),... %R2.*1000 for um->nm
    'XData',[0 size(G, 2).*cvFt],...
    'YData',[0 size(G, 1).*cvFt])
hold all
%% Vol3d plot - display tweaks
% cmap = colormap(gca);
% amap = alphamap(gca);
%load('Z:\RicardoMongeN\IX-73 collections\2023_04_20\F1_a1_zMap_2kfpSt_5uL_10mW_30ms_400g\cmap.mat');
%load('Z:\RicardoMongeN\IX-73 collections\2023_04_20\F1_a1_zMap_2kfpSt_5uL_10mW_30ms_400g\amap.mat');
load('Z:\RicardoMongeN\IX-73 collections\2023_04_20\F1_a1_zMap_2kfpSt_5uL_10mW_30ms_400g\cmap.mat');
colormap(cmap); %brighten(0.3);%shading interp;
axis image;colorbar;%('East','Color','w');%colormap hot;
set(gca,'Color','w')
set(gca,'GridColor','k')
%imcontrast
%alphamap([0 linspace(0.1, 0, 255)]);


%%
view(45,45)
title(' ')%strcat('CSP Localizations Map'), 'FontSize',16)
ylabel('Y (nm)','FontSize',20, 'Rotation',40) %use 40
xlabel('X (nm)','FontSize',20, 'Rotation',-40) %use -40
zlabel('Z (nm)','FontSize',20)
set(gca,'LineWidth',2);

ax = gca;
ax.FontSize = 18;

% set(gca,'Color','k')
% set(gca,'GridColor','w')
hold on
grid on

%% Vold3d replotting
view(45,45)
% zlim([-500 6000])
% xlim([500 7000])
% ylim([1000 7500])

zlim([-500 2500])
xlim([850 5850])
ylim([1500 6500])
%%

xticklabels("");yticklabels("");zticklabels("");
xlabel("");ylabel("");zlabel("");

colorbar off