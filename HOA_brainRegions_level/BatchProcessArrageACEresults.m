%% batch process ACE results examination
clear all;
path.data='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info';

% Stats={'AAL','BA','BNA'};
% filelist=spm_select('List',path.data,'twins_data_CBF_.*forACE_UnivAE.csv');
filelist='twins_data_CBF_HOA_whole_withSmooth_forACE_UnivAE.csv';
filelist=cellstr(filelist);

for i=1:length(filelist)
    
    AE_dataDir=fullfile(path.data,filelist{i}); %twins_data_CBF_AAL_withoutSmooth_forACE_UnivAE.csv
    ICC_dataDir=strrep(filelist{i},'twins_data_CBF','icc_twin');   % icc_twin_corr_BA_withSmooth_forACE.csv
    ICC_dataDir=strrep(ICC_dataDir,'_UnivAE','');
    ICC_dataDir=fullfile(path.data,ICC_dataDir);
%     ICC_dataDir=fullfile(path.data,['icc_CBF_' Stats{i} '.csv']);
   temp_name=strrep(filelist{i},'twins_data_CBF_','');
    temp_name=strrep(temp_name,'_forACE_UnivAE.csv','');
    
    AIC_pertOutputDir=fullfile(path.data,['ResultsArr_' temp_name '.mat']);
    OutputDir=fullfile(path.data,['ResultsArr_sigTest_' temp_name '.csv']);
ArrageACEresults(AE_dataDir,ICC_dataDir,AIC_pertOutputDir,OutputDir)
end