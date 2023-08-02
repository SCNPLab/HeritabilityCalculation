%% arrange data
clear all;
% path.data='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\normalizedCBF\CBFtwins';  % without smooth
path.data='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\normalizedCBF\smoothedCBFtwins';  % with smooth
path.output='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info';
path.roi='F:\IPCAS_TWIN\CBF\info\Atlases\Reslice_Atlases';

% atlas_set={'AAL','BNA','BA','HOAc'};
atlas_set={'AAL','BNA','BA','HOA_whole' ,'HOAc'};

demog_data=readtable('F:\IPCAS_TWIN\CBF\ASL_template.xlsx','Sheet','AAL');
demog_data=demog_data(:,1:5);

ROIset=spm_select('FPList',path.roi,'.nii');
ROIset=cellstr(ROIset);

filelist=spm_select('FPList',path.data,'.nii');
filelist=cellstr(filelist);


for ROIord=4%1:length(ROIset)
    
    RefImage=load_nii(ROIset{ROIord});
    NumRegions=unique(RefImage.img);
    NumRegions(NumRegions==0)=[];
    
    CBF.(atlas_set{ROIord})=zeros(length(filelist),length(NumRegions));
    for i=1:length(filelist)
    
        for k=1:length(NumRegions)
            temp_image=load_nii(filelist{i});
            temp_mean=mean(temp_image.img(RefImage.img==NumRegions(k)),'all');
            
            CBF.(atlas_set{ROIord})(i,k)=temp_mean;
        end
    end
    
    temp_data=CBF.(atlas_set{ROIord});
    temp_data=array2table(temp_data); 
        
     VariableNum=size(temp_data,2);
variableNames=cell(1,VariableNum);

for volOrd=1:VariableNum
    variableNames{volOrd}=['CBF_' num2str(volOrd)];
end

    temp_data.Properties.VariableNames=variableNames;
    outputT=[demog_data,temp_data];
    % write output
%     writetable(outputT,fullfile(path.output,['CBF_' atlas_set{ROIord} '_withoutSmooth_forACE.csv']),'Delimiter',',');   
    writetable(outputT,fullfile(path.output,['CBF_' atlas_set{ROIord} '_withSmooth_forACE.csv']),'Delimiter',',');   
    fprintf('\n %s finished!',atlas_set{ROIord});
        
end       
        
%   save(fullfile(path.output,'CBFtwins_withSmooth.mat'),'CBF');      
% % for i=1:length(atlas_set)
% %     
% % %     [temp_data,~,~]=xlsread(fullfile(path.data,'ROIcbf_nsr.xlsx'),atlas_set{i});
% % %     temp_data=array2table(temp_data);
% %  temp_data=load(fullfile(path.data,['ROISignals_ROISignal_' atlas_set{i} '_CBFtwins.mat']));
% %  temp_data=temp_data.ROISignals;
% %  temp_data=array2table(temp_data);
% % %     temp_data=array2table(temp_data);
% %         % variable names
% % % [temp,~,~]=xlsread(fullfile(path.data,'ROIcbf_nsr.xlsx'),'AAL');
% % VariableNum=size(temp_data,2);
% % variableNames=cell(1,VariableNum);
% % 
% % for k=1:VariableNum
% %     variableNames{k}=['CBF_' num2str(k)];
% % end
% % 
% %     temp_data.Properties.VariableNames=variableNames;
% %     outputT=[demog_data,temp_data];
% %     writetable(outputT,fullfile(path.output,['CBF_' atlas_set{i} '_withSmooth_forACE.csv']),'Delimiter',',');
% %     
% % 
% % 
% % end
% %     
    
