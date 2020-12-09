%% LRG_SuperRes_Read_Data
% Read data from .sif files and save it as a .mat file for use with the
% analysis code
function LRG_SuperRes_Read_Data(e)
cd(e.path)
filename=strcat(e.filename,'.sif');
Data1=AndorReadAllFrames(filename,e.path,e.startframe,e.stopframe);
save(strcat(e.filename,'.mat'),'Data1')
return