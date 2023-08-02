%% plot heritability of CBF Dang 20211103
clear all;close all;
path.data='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info\HOA_20221006';
path.output='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info\HOA_20221006';
path.atlas='F:\IPCAS_TWIN\CBF\info\Atlases\Reslice_Atlases';
path.label='F:\IPCAS_TWIN\CBF\info\Info';
% wholeGroupSet={'AAL','BA','BNA'};
% fileset=spm_select('List',path.data,'twins_data_CBF.*forACE_UnivAE.csv');
% fileset='twins_data_CBF_HOA_whole_withSmooth_forACE_UnivAE.csv';
fileset='BestModelInfo_3.mat';
fileset=cellstr(fileset);

% label AAL
region_labels.aal=load(fullfile(path.label,'aal_Labels.mat'));
region_labels.aal=region_labels.aal.Reference(:,1);
region_labels.aal(1,:)=[];
% label BA
region_labels.ba=load(fullfile(path.label,'Brodmann_YCG_Labels.mat'));
region_labels.ba=region_labels.ba.Reference(:,1);
region_labels.ba(1,:)=[];
% label HOVc
region_labels.hov=load(fullfile(path.label,'HarvardOxford-cort-maxprob-thr25-2mm_YCG_Labels.mat'));
region_labels.hov=region_labels.hov.Reference(:,1);
region_labels.hov(1,:)=[];

region_labels.hov_sub=load(fullfile(path.label,'HarvardOxford-sub-maxprob-thr25-2mm_YCG_Labels.mat'));
region_labels.hov_sub=region_labels.hov_sub.Reference(:,1);
region_labels.hov_sub(1,:)=[];

region_labels.hov_whole=[region_labels.hov;region_labels.hov_sub];

% label BNA
region_labels.bna=readtable(fullfile(path.label,'Yeo_BN.xlsx'));
region_labels.bna=table2cell(region_labels.bna(:,3));

ROIlabels={'aal','ba','bna','hov'};


% ROIlist=spm_select('List',path.atlas,'.nii');
% ROIlist=cellstr(ROIlist);
% ROIset={ROIlist{1};ROIlist{3};ROIlist{2};ROIlist{4}};
  ROIset={'HOC_cortical_subcortical.nii'};
 for wholeGroupSetOrd=1:length(fileset)
   
%     temp_name=fileset{wholeGroupSetOrd};
%     temp_name=strrep(temp_name,'twins_data_CBF_','');
%     temp_name=strrep(temp_name,'_forACE_UnivAE.csv','');
    temp_name='HOA_whole';
%     
%     heritaData=readtable(fullfile(path.data,fileset{wholeGroupSetOrd}));
%     SigInfo=readtable(fullfile(path.data,['ResultsArr_sigTest_' temp_name '.csv']));
      
heritaData_raw=load(fileset{wholeGroupSetOrd});
heritaData_raw=heritaData_raw.output;



%     SigInfo=SigInfo.PPMaic;
    SigInfo=heritaData_raw.best_model_para(:,13);
    SigInfo_C=heritaData_raw.best_model_para(:,14);
    
%     herita_a=heritaData.a;
%     herita_e=heritaData.e;
    herita_a=heritaData_raw.best_model_para(:,1);
    herita_c=heritaData_raw.best_model_para(:,4);
    
%     halfNum=length(herita_a)/2;
%     herita_a(1:halfNum)=[];
%     herita_e(1:halfNum)=[];
    
%      combinedA_E=[herita_a,herita_e];
    combinedA_E=herita_a;
    
    Sig_status=SigInfo>0.9;
    
    Sig_status_C=SigInfo_C>0.9;
    %% for create a table
    select_best_model_para=heritaData_raw.best_model_para(Sig_status,:);
    select_best_model=heritaData_raw.modelNames(Sig_status,:);
    select_region_labels=region_labels.hov_whole(Sig_status,:);
%     Sig_status=double(Sig_status);
    
    herita_a_aft_ppm=herita_a.*Sig_status;
    herita_c_aft_ppm=herita_c.*Sig_status_C;
%    herita_a_aft_ppm=herita_a(Sig_status);
% % %     % PPM
% % %     PPM=SigInfo;
% %     
% % %     PPM(1:halfNum)=[];
% %     
% % %     LabelsSig=cell(1,length(PPM));
% % %     for k=1:length(PPM)
% % %         if PPM(k)>0.9
% % %         LabelsSig{1,k}='*';
% % %         else
% % %         LabelsSig{1,k}='';
% % %         end
% % %     end
% % %     LabelsSig=string(LabelsSig);
% %     LabelsSig=cell(1,length(herita_a_aft_ppm));
% %     for k=1:length(herita_a_aft_ppm)
% %         if herita_a_aft_ppm(k)>0
% %         LabelsSig{1,k}='*';
% %         else
% %         LabelsSig{1,k}='';
% %         end
% %     end
% %     LabelsSig=string(LabelsSig);
% %     % plot figure
% %      figure(1)
% % %     b=bar(combinedA_E,'stacked');
% % %     b=bar(combinedA_E);
% %     b=bar(herita_a_aft_ppm);
% %     xtips1 = b(1).XEndPoints;
% %     ytips1=b(1).YEndPoints;
% %      for yOrd=1:length(ytips1)
% %         if ytips1(yOrd)<0
% %             ytips1(yOrd)=ytips1(yOrd)-0.2;
% %         end
% %     end
% %     text(xtips1,ytips1,LabelsSig,'HorizontalAlignment','center',...
% %     'VerticalAlignment','bottom');
% % %     xticks([1:48]);
% % %     xtickangle(45);
% % %     xticklabels(DTInames);
% %         remOrd=ceil(wholeGroupSetOrd/2);
% %          which_labels=region_labels.('hov_whole');
% %          which_labels=which_labels(Sig_status);
% %          for labOrd=1:length(which_labels)
% %              which_labels{labOrd}=strrep(which_labels{labOrd},'_','-');
% %          end
% %          
% %          xticks(1:length(herita_a_aft_ppm));
% %          xtickangle(45);
% %          xticklabels(which_labels);
% %     set(gcf,'WindowState','Maximized');
% %     set(gca,'FontSize',12);
% % %     title(temp_name);
% %      saveas(1,fullfile(path.output,['SigHerita_aftPPM_' temp_name '.png']),'png');
% %       close all;
% %     
      
%       RefImageDir=fullfile(path.atlas,ROIset{1});
     % write heritability info into nii
% %        if wholeGroupSetOrd<=2
                RefImageDir=fullfile(path.atlas,ROIset{1});
              
% %        end
             clear V 
     V=spm_vol(RefImageDir);
            RefImage=spm_read_vols(V);
            NumRegions=unique(RefImage);
            NumRegions(NumRegions==0)=[];
             for i=1:length(NumRegions)
%                 RefImage(RefImage==NumRegions(i))=herita_a(i); % herita_a_aft_ppm
                RefImage(RefImage==NumRegions(i))=herita_a_aft_ppm(i);
%              RefImage(RefImage==NumRegions(i))=herita_c_aft_ppm(i);
             end
            V.fname=fullfile(path.output,[temp_name '_CBF_heritability_aftPPM_3.nii']);
%             V.fname=fullfile(path.output,[temp_name '_CBF_sharedEnv_aftPPM.nii']);
            V.dt(1)=16;
%             spm_write_vol(V,RefImage);
%              fprintf('\n %s finished',temp_name);
end
    
    
    