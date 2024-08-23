%%  Shifts in XY and Z for Bead Scans / Grouped Sites
%7/11/22
clear
close all
clc

%Uses only LRG_SuperRES output data, for Z-scans

% load('Z:\RicardoMongeN\IX-73 collections\2021_05_24\b1, 532, z scan, 13.6mW, 100xg, 100ms\localizations_filtered_2021_08_11\Output.mat')
startloc='Z:\RicardoMongeN\IX-73 collections\2022_10_27\';
numFiles =1; %how many times to loop over for single files
multiSelect = 'false';%allow for multiselect - use if all in single folder

fileNames = "";
filePaths = "";

if strcmp(multiSelect,'true')==0
    for i=1:numFiles
        [file, path] = uigetfile(startloc,...
            'Select a analyzed MAT image file')
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

%% Run XY Tracks
close all

StDev = struct([]);

SHX = []; %sigma X from mean, for all
SHY = [];%sigma Y from mean, for all
FWHMx = []; %full width half-max of fitted Gaussians, for all
FWHMy = [];

avGWx = [];
avGWy = [];

for j=1:numFiles

%% Get Z Positions OR Specify if Stationary
% 
addpath(filePaths(j+1));
% 
% Digits = 3; %significant digits in Z measurements
% stationary = 1;% 0=No,1=Yes
% 
% 
% fileID=fopen(filePaths(j+1),'rt'); %open the file
% 
% if strcmp(multiSelect,'true')==1
%     if j==1
%         [file,path] = uigetfile(startloc,...
%        'Select a Metadata file','.txt');
%         txtFile =  file;
%     end
% else
%     
%     txtFile = erase(fileNames(j+1),'.ome.tif_analyzed.mat');
%     txtFile = strcat(txtFile,'_metadata.txt');
% end
% 
% A = fileread(txtFile);
% 
% startPat ='"ZPositionUm": '; %mark Start pattern to search for in File array
% endPat = ','; %mark End pattern to search for in File array
% 
% newA = extractBetween(A,startPat,endPat); %generate cell array, containing strings located between the Start and End patterns.
% B = newA'; %take transpose of cell array
% 
% C = [Inf 1]; %initialize new array containing numerical form of B
% for i =1:length(B)
%     D = cell2mat(B(i)); %convert the cell arrays to ordinary arrays
%     C(i) = str2num(D);  %convert str arrays to numerical values
%     if stationary == 1
%         C(i) = 0;
%     end
% end
% zList = round(C',Digits,'significant'); %round position array to the correct number of significant digits
% % zList = zList.*1000; %convert to nm
% 
% clear A B newA C D

%% Load and Organize Output Data
clc

load(fileNames(j+1));


pixelSize = e.pixelSize;%nm

ogLocatStore = LocatStore;
frameList =1:1:size(ogLocatStore,2); %list of unique frame numbers
nFrames = length(frameList); %how many frames were particles ID'd in  


%% GroupLocat Analyzis
clc;%close all


for i=1:size(GroupLocat,2)
    sigX = mean(GroupLocat(i).RawSites(:,2)-GroupLocat(i).Centroid(2))*e.pixelSize;%std(GroupLocat(i).RawSites(:,2))*pixelSize;%in nm
    sigY = mean(GroupLocat(i).RawSites(:,1)-GroupLocat(i).Centroid(1))*e.pixelSize;%std(GroupLocat(i).RawSites(:,1))*pixelSize;%in nm
    sigXY = sqrt(sigX^2+sigY^2);
    
    StDev(j).SigmaX(i,:) = sigX;
    StDev(j).SigmaY(i,:) = sigY;
    StDev(j).SigmaXY(i,:) = sigXY;
    StDev(j).Int(i,:) = mean(GroupLocat(i).RawSites(:,4));
    StDev(j).SigmaI(i,:) = std(GroupLocat(i).RawSites(:,4));  
    
    
    figure(1);hold on; plot(StDev(j).Int(i,:),StDev(j).SigmaX(i,:),'*');
    xlabel('Intensity (AU)');ylabel('Sigma X (nm)');
    figure(2);hold on; plot(StDev(j).Int(i,:),StDev(j).SigmaY(i,:),'*');
    xlabel('Intensity (AU)');ylabel('Sigma y (nm)');
    
    SHX = [SHX, transpose(abs((GroupLocat(i).RawSites(:,2))-mean(GroupLocat(i).RawSites(:,2)))*e.pixelSize)];
    SHY = [SHY,  transpose(abs((GroupLocat(i).RawSites(:,1))-mean(GroupLocat(i).RawSites(:,1)))*e.pixelSize)];
end

tempX = [];
tempY = [];
for i=1:size(LocatStore,2)
    FWHMx = [FWHMx, transpose(LocatStore(i).PSFfinal(:,4).*e.pixelSize)];
    FWHMy = [FWHMy, transpose(LocatStore(i).PSFfinal(:,3).*e.pixelSize)];
    tempX = [tempX, mean(transpose(LocatStore(i).PSFfinal(:,4).*e.pixelSize))];
    tempY = [tempY, mean(transpose(LocatStore(i).PSFfinal(:,3).*e.pixelSize))];
end
cT = 0;
cF = 100;
for i=1:size(tempX,2)/cF;
    avGWx = [avGWx, mean(tempX(cT*cF+1:i*cF))];
    avGWy = [avGWy, mean(tempY(cT*cF+1:i*cF))];
    cT = cT+1;
end

end
% FWHMx = FWHMx(FWHMx(:)< 240);
% FWHMy = FWHMy(FWHMy(:)< 240);

figure(3);hold on; histogram(SHX);xlim([0, 200]);
xlabel('X Deviation from Mean (nm)');ylabel('Beads over all frames');
figure(4);hold on; histogram(SHY);xlim([0, 200]);
xlabel('Y Deviation from Mean (nm)');ylabel('Beads over all frames');
figure(5);hold on; histogram(FWHMx);
xlabel('X Gaussian Width (nm)');ylabel('Beads');
figure(6);hold on; histogram(FWHMy);
xlabel('Y Gaussian Width (nm)');ylabel('Beads');

%%
figure(7);hold on; plot(avGWx);
xlabel('Frame'); ylabel('X Mean PSF Width (nm)');
figure(8);hold on; plot(avGWy);
xlabel('Frame'); ylabel('Y Mean PSF Width (nm)');


%% Save current File Output
save(strcat(startloc,'XYtracks.mat'),'StDev')

%% Open XYtracks for multiple files
clear;close all;clc
startloc='Z:\RicardoMongeN\IX-73 collections\2022_07_25\High\0um\';
numFiles =1;

fileNames = "";
filePaths = "";
for i=1:numFiles
    [file, path] = uigetfile(startloc,...
        'Select a XYtracks MAT image file')
    fileNames = [fileNames, file];
    filePaths = [filePaths, path];
end
disp('All Files Selected')

%% Extract sigmaX, sigmaY, I, sigmaI for all files

sigmaX =[];
sigmaY = [];
sigmaXY = [];
sigmaI = [];
meanI = [];

for i=1:numFiles
   load(strcat(filePaths(i+1),fileNames(i+1)));
   for j=1:size(StDev,2)
       sigmaX = [sigmaX, transpose(StDev(j).SigmaX(:))];
       sigmaY = [sigmaY, transpose(StDev(j).SigmaY(:))];
       sigmaXY = [sigmaXY, transpose(StDev(j).SigmaXY(:))];
       sigmaI = [sigmaI, transpose(StDev(j).SigmaI(:))];
       meanI = [meanI, transpose(StDev(j).Int(:))];
   end
    
end

resX = [mean(sigmaX),std(sigmaX)];
resY = [mean(sigmaY),std(sigmaY)];
resXY = [mean(sigmaXY),std(sigmaXY)];

disp(strcat("SigmaX = ",string(resX(1)),"+/- ",string(resX(2))," nm"))
disp(strcat("SigmaY = ",string(resY(1)),"+/- ",string(resY(2))," nm"))
disp(strcat("SigmaXY = ",string(resXY(1)),"+/- ",string(resXY(2))," nm"))
%% Plot and Fit 1/sqrt(N)
close all;%clc
%figure(1);hold on;plot(0:1:7000,1000./sqrt(0:1:7000),'-')
warning off
%Sigma X Fits
    y = sigmaX;
    x = meanI;
    xErr = sigmaI;
    
    f = @(b,x) b(1)./(sqrt(x)) + b(2); %model
    beta0 = [1000;1];
    
    OLS = @(b) sum((f(b,x) - y).^2); % Ordinary Least  Squares cost function
    B1 = fminsearch(OLS, beta0);  % Estimate Parameters
    Rsq2 = 1 - sum((y - f(B1,x)).^2) / sum((y - mean(y)).^2);

%         opts = statset('MaxIter',200);
%         mdl = fitnlm(x,y,f,beta0,'Options',opts);
%         B3 = mdl.Coefficients.Estimate;
%         Rsq3 = mdl.Rsquared.Adjusted;

    figure();hold on;
    plot(x, y, 'or')
    Xrange = linspace(min(x),max(x));
    plot(Xrange, f(B1,Xrange), '-')
    xlabel('Intensity (AU)')
    ylabel('X Sigma (nm)')
    
    text(mean(x), 1.5*mean(y),...
    sprintf('f(x) = %.1f/sqrt(x) + (%.1f)%', B1),...
    'Color','b')
   
%Sigma Y Fits
    y2 = sigmaY;
   
    f2 = @(b,x) b(1)./(sqrt(x)) + b(2); %model
    beta0 = [1000;1];
    
    OLS = @(b) sum((f2(b,x) - y2).^2); % Ordinary Least Squares cost function
    B2 = fminsearch(OLS, beta0);  % Estimate Parameters
    Rsq2 = 1 - sum((y2 - f2(B2,x)).^2) / sum((y2 - mean(y2)).^2);

%         opts = statset('MaxIter',200);
%         mdl = fitnlm(x,y,f,beta0,'Options',opts);
%         B3 = mdl.Coefficients.Estimate;
%         Rsq3 = mdl.Rsquared.Adjusted;

    figure();hold on;
    plot(x, y2, 'or')
    Xrange = linspace(min(x),max(x));
    plot(Xrange, f(B2,Xrange), '-')
    xlabel('Intensity (AU)')
    ylabel('Y Sigma (nm)')
    
    text(mean(x), 1.5*mean(y2),...
    sprintf('f(x) = %.1f/sqrt(x) + (%.1f)%', B2),...
    'Color','b')
    
    
%Sigma XY Fits
    y3 = sigmaXY;
   
    f3 = @(b,x) b(1)./(sqrt(x)+b(2)) +b(3); %model
    beta0 = [1000;1;1];
    
    OLS = @(b) sum((f3(b,x) - y3).^2); % Ordinary Least Squares cost function
    B3 = fminsearch(OLS, beta0);  % Estimate Parameters
    Rsq2 = 1 - sum((y3 - f3(B3,x)).^2) / sum((y3 - mean(y3)).^2);

%         opts = statset('MaxIter',200);
%         mdl = fitnlm(x,y,f,beta0,'Options',opts);
%         B3 = mdl.Coefficients.Estimate;
%         Rsq3 = mdl.Rsquared.Adjusted;

    figure();hold on;
    plot(x, y3, 'or')
    Xrange = linspace(min(x),max(x));
    plot(Xrange, f3(B3,Xrange), '-')
    xlabel('Intensity (AU)')
    ylabel('XY Sigma (nm)')
    
    text(mean(x), 1.5*mean(y3),...
    sprintf('f(x) = %.1f/(sqrt(x)+(%.1f)) + (%.1f)%', B3),...
    'Color','b')
