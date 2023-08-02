%% arrange data
clear all;
% path.data='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\normalizedCBF\CBFtwins';  % without smooth
path.data='H:\IPCAS_TWIN\CBF\SmoothedData\20220620\normalizedCBF\smoothedCBFtwins';  % with smooth
path.output='H:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info\average4Cortex_HOV\extractedCBF_divBilateral';
if ~exist(path.output,'dir')
    mkdir(path.output);
end
path.roi='H:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info\average4Cortex_HOV\reference_images\DivBilateral';

% atlas_set={'AAL','BNA','BA','HOAc'};
% atlas_set={'AAL','BNA','BA','HOA_whole' ,'HOAc'};
% atlas_set=spm_select('List',path.roi,'.nii');
% atlas_set=cellstr(atlas_set);


demog_data=readtable('H:\IPCAS_TWIN\CBF\ASL_template.xlsx','Sheet','AAL');
demog_data=demog_data(:,1:5);

ROIset=spm_select('List',path.roi,'.nii');
ROIset=cellstr(ROIset);
ROInames=cell(length(ROIset),1);

filelist=spm_select('FPList',path.data,'.nii');
filelist=cellstr(filelist);

CBF_output=zeros(length(filelist),length(ROIset));


for ROIord=1:length(ROIset)
    
%     RefImage=load_nii(ROIset{ROIord});
%     NumRegions=unique(RefImage.img);
  RefImage_V=spm_vol(fullfile(path.roi,ROIset{ROIord}));
   RefImage_img=spm_read_vols(RefImage_V);
%     NumRegions=unique(RefImage.img);
%     NumRegions(NumRegions==0)=[];
  ROIname=ROIset{ROIord};
  ROIname=strrep(ROIname,'.nii','');
  ROIname=strrep(ROIname,'-','_');
  ROInames{ROIord}=ROIname;

  
    for i=1:length(filelist)
    
      
            temp_image_V=spm_vol(filelist{i});
            temp_image=spm_read_vols(temp_image_V);
            temp_mean=mean(temp_image(RefImage_img>0.2),'all');
            
           CBF_output(i,ROIord)=temp_mean;
      
    end
    
   
        fprintf('\n ROI %s finished',ROIname);
end       
        
   CBF_output_mat.mat=CBF_output;
   CBF_output_mat.ROInames=ROInames;
   save(fullfile(path.output,'CBF.mat'),'CBF_output_mat');
   
   % save as csv for ACE model
    CBF_output=array2table(CBF_output); 
        
     VariableNum=size(CBF_output,2);
variableNames=cell(1,VariableNum);

for volOrd=1:VariableNum
%     variableNames{volOrd}=['CBF_' num2str(volOrd)];
   variableNames{volOrd}=ROInames{volOrd}; 
%    variableNames{volOrd}=strrep(variableNames{volOrd},'-','_');
end

    CBF_output.Properties.VariableNames=variableNames;
    outputT=[demog_data,CBF_output];
    % write output
%     writetable(outputT,fullfile(path.output,['CBF_' atlas_set{ROIord} '_withoutSmooth_forACE.csv']),'Delimiter',',');   
    writetable(outputT,fullfile(path.output,'CBF_cortex_arterySubregions_withSmooth_forACE.csv'),'Delimiter',',');   
  
    
