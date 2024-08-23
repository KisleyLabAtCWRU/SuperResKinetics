%% Sphere Test_V2
%3D volume plot of Super-Res zScan data

%Plot localizations for a long z-scan
%First, obtain Log of Z positions
%RicardoMN - KL
clear;close all;clc;

[file,path] = uigetfile('*.txt');
addpath(path);

%StepSize = 0.1; %in um  
Digits = 3; %significant digits in measurements

fileID=fopen(path,'rt'); %open the file
 
A = fileread(file);

startPat ='"ZPositionUm": '; %mark Start pattern to search for in File array
endPat = ','; %mark End pattern to search for in File array

newA = extractBetween(A,startPat,endPat); %generate cell array, containing strings located between the Start and End patterns.
B = newA'; %take transpose of cell array

C = [Inf 1]; %initialize new array conta ining numerical form of B
for i =1:length(B)
    D = cell2mat(B(i)); %convert the cell arrays to ordinary arrays
    C(i) = str2num(D);  %convert str arrays to numerical values
end
R= round(C,Digits,'significant'); %round position array to the correct number of significant digits
% R = R.*1000; %convert to nm

clearvars -except R
%% Plot Localizations 
%RMN 3/11/2022

%Use LRG_SuperRes to identify particles over many
%frames (bining) - Or use and skip to 3D plotting (bottom)
%***Load ZLog file (made with "MakeZLog" script) AND 3D Localization Outputs for 3D Plot Comparisons
close all


aList = [1]; %corresponds to files used

list2Use = aList; %tell what list to use as reference for the number of files needed
numFiles = length(list2Use);


fileNames = "";
filePaths = "";
for i=1:numFiles
    [file, path] = uigetfile('.tif','Select a TIF image file')
    fileNames = [fileNames, file];
    filePaths = [filePaths, path];
end



%% User Defined Parameters - Run section for 2D and 3D
% General Parameters
e.runsuperres='true'; % (true or false) false = find kinetics using diffraction limited data. true = find kinetics using superlocalized data
e.xmin=1;
e.xmax=512;
e.ymin=1;
e.ymax=512;
e.pixelSize=160; % nm

% Data Parameters
e.loadsif='true'; % true = load data from a .sif file (if this is the first analysis for instance) false = data will be loaded from a .mat file

% Particle Identification Parameters
e.BackgroundThreshold = 1; %multiplicative minimum treshold over local background for particle ID (before adding sigma)
e.selectROI = 'false'; %prompts uuser to draw an ROI for each file
e.wasCut = 'false';

wasDone = 0;% if Sphere_v2 analysis was done before = 1, otherwise =0;

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
e.FinalLocatSigma=0.5;%loc accuracy-the standard deviation of the single binding site - used for model 2D peak (0.5)
e.nevent=3;%events-minimum peak intensity of the spot - to distinguish specific and non-specific intrs (5 events)
e.SRCorrFactor=0.6; %shape-the minimum cross correlation factor between a spot and the standard Gaussian peak (0.6)
e.FinalLocatThresh=3;%search radius - how many FinalLocatSigma are used to rule out spfc vs non-spfc (3 sigma)

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

Analyzed = struct([]);
DwellSt = {};

DwellTable = struct([]);

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
    %e.startframe = 1;
    e.stopframe=fileInfo(3);
    e.xmax=fileInfo(1);
    e.ymax=fileInfo(2);
    
    disp('Starting')
    
    disp('Generating and Cutting Data as .mat as needed')  
    
    %find unique indeces of R(Z position vector)
    [R2,ia,ic] = unique(R(1,1:fileInfo(3))); %[unique positions, R2=R(ia),R=R2(ic)]
    
    %Process and plot data for specific ranges
    
    for j=1:(size(R2,2)-1); %j is unique Z position in list R2
        Analyzed(j).Data = [[],[]];
        
        
        %Dynamic frame ranges
        e.startframe = ia(j);
        e.stopframe = ia(j+1)-1;%frame right before where position changed
        e.nframes=e.stopframe-e.startframe+1;
        %e.kinetiC = 'true';
        %e.loadsif = 'false';
        
        if wasDone == 0;
            e.loadsif = 'true';
            disp('Analyzing current slice')
            RM_LRG_SuperRes_Run(e);  %comment out if not running localization code,
            %and simply load in the files
            disp('Loading Current Analyzed slice')
            load(strcat(e.filename,'_analyzed.mat'));
            if strcmp(e.kinetiC,'true')==1;
                save(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'),'LocatStore','e','Data1','GroupLocat','Ensemble'); 
            else 
                save(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'),'LocatStore','e','Data1'); 
            end
            load(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'));
            
        else
            disp(strcat("Loading ",num2str(R2(j)),'um Analyzed File'))
            
            load(strcat(e.filename,"_",num2str(R2(j)),'um_analyzed.mat'));
            e.loadsif = 'false';
        end
  
        Fm = e.startframe:1:e.stopframe; %frame numbers to match limits
        
        %Making XxYxZ array for vol3d function
        Analyzed(j).Data(:,:) = LRG_SuperRes_GenerateSR(LocatStore,e);
        Analyzed(j).Sites = [[],[]];
        Analyzed(j).siteDwell = [[],[]];
        Analyzed(j).siteLocs = [[],[]];
       
        if strcmp(e.kinetiC,'true')==1;     
        for k = 1:size(GroupLocat,2)
           Analyzed(j).Sites(k,:) = GroupLocat(k).Centroid(:);%site [y,x,dx,dy] per Z
%            Analyzed(j).siteDwell(k,:) = GroupLocat(k).Dwell(:); %site dwell times per Z
           Analyzed(j).siteLocs(k,:) = size(GroupLocat(k).RawSites,1); %number of locs per site at Z        
           
        end
        end
        %Store x and y Gaussian Widths for each particle at current R(z)
        %value
%         Analyzed(j).xGaussWidth = zeros(size(LocatStore,2),1);
%         Analyzed(j).yGaussWidth = zeros(size(LocatStore,2),1);
        kCount = 1; %particle count for current Z position R2(j)
        for k = 1:size(LocatStore,2) %frame k, in file i, for j Z position, for molecule h
            if isempty(LocatStore(k).PSFfinal)==0
                for h =1:size(LocatStore(k).PSFfinal,1)
                    Analyzed(j).xGaussWidth(kCount) = LocatStore(k).PSFfinal(h,4);
                    Analyzed(j).yGaussWidth(kCount) = LocatStore(k).PSFfinal(h,3);
                    Analyzed(j).Intensity(kCount) = LocatStore(k).PSFfinal(h,5);             
                    kCount = kCount+1;
                end
            end
        end
        
        %Store [mean,stdev] Gauss Width and Intensity at each R(z),...
        % jth position
        Analyzed(j).mX = [mean(Analyzed(j).xGaussWidth(:)),std(Analyzed(j).xGaussWidth(:))];
        Analyzed(j).mY = [mean(Analyzed(j).xGaussWidth(:)),std(Analyzed(j).xGaussWidth(:))];
        Analyzed(j).mI = [mean(Analyzed(j).Intensity(:)),std(Analyzed(j).Intensity(:))];
        
        if strcmp(e.kinetiC,'true')==1;
        Analyzed(j).enDwell = Ensemble.Dwell; % Dwell for each Z slice
        
        
        %cumuldist for ensemble dwell at each Z slice
        if isempty(Ensemble.Dwell)==0
            [Tm, Pt] = cumuldist(Ensemble.Dwell, unique(Ensemble.Dwell));

            DwellTable(j).Pt = Pt;
            DwellTable(j).Tm = Tm;
         
            DwellSt{2*j,1} = Pt;%time
            DwellSt{2*j+1,1} = Tm;%prob
        else
            DwellSt{2*j,1} = [];
            DwellSt{2*j+1,1} = [];
        end
        
        end
        %Localization Distribution Arrays
        numLocs = [];
        for h =1:size(Fm,2)
            numLocs = [numLocs, size(LocatStore(h).PSFfinal,1)];
        end
        Analyzed(j).NumLocs = sum(numLocs);
        Analyzed(j).LocsPerFrame = sum(numLocs)/size(Fm,2);
      
    end
    
    
    %writecell(DwellSt,strcat((e.path),'DwellSt.csv'));
    disp('Done')

end

%writecell(DwellSt,strcat((e.path),'DwellSt.csv'));



%%
%cmap = colormap(gca);
%amap = alphamap(gca);
% load('Z:\RicardoMongeN\Chiral - Working Outline\ImageData\3D\cmap.mat');
% load('Z:\RicardoMongeN\Chiral - Working Outline\ImageData\3D\amap.mat');
%colormap(cmap); %brighten(0.3);%shading interp;
%axis image;colorbar;%('East','Color','w');%colormap hot;
%imcontrast
%alphamap([0 linspace(0.1, 0, 255)]);

%% Test vol3d
%addpath('Z:\RicardoMongeN\MATLAB')
%clearvars -except Analyzed R2 e file fileInfo fileNames filePaths grad rEx
G = [[],[],[],[]];
%grad2 =  colorGradient([0,0,0],[0.2,0.8,0.3],size(R2,2));
%load mri.mat;
disp('Arranging SR Image Data')
for i=1:floor((size(R2,2)-1));
    G(:,:,1,i)=Analyzed(i).Data(:,:); 
end
disp('Done')
%hold on;imagesc(sum(squeeze(G),3));colormap(cmap);alphamap(amap);
% set(gca,'Color','w')
% set(gca,'GridColor','k')
save(strcat(filePaths(2),fileNames(2),"_3DmapData.mat"),'G','R2','e');
disp('Saved 3D Map Data');
%%  Vold3d plotting
clear Analyzed
addpath('Z:\RicardoMongeN\MATLAB')
figure()

G = squeeze(G());
%G = G(:,:,1:5);
%clim([0 25])
vol3d('Cdata',squeeze(G()),...
    'ZData',R2.*1000,... %R2.*1000 for um->nm
    'XData',[0 size(G, 2).*(e.pixelSize/e.nzoom)],...
    'YData',[0 size(G, 1).*(e.pixelSize/e.nzoom)])
hold all
%% Overlay Sorption Sites Located at Each Z_start
% for i=1:size(Analyzed,2)
%     for j=1:size(Analyzed(i).Sites,1)
%         scatter3(Analyzed(i).Sites(j,2).*(e.pixelSize)...
%             ,Analyzed(i).Sites(j,1).*(e.pixelSize)...
%             ,R2(i).*1000,'r')
%     end
% end

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
% view(45,45)
% % zlim([-500 6000])
% % xlim([500 7000])
% % ylim([1000 7500])
% 
% zlim([-500 6500])
% xlim([5500 10000])
% ylim([5500 10000])
%%

xticklabels("");yticklabels("");zticklabels("");
xlabel("");ylabel("");zlabel("");

colorbar off
