%% display the best model 
%% batch process ACE results examination
clear all;
% path.data='F:\IPCAS_TWIN\dFC\dFC\Results\ProbabilityOfCoactivation\ResultsOfDifferentCombGammaOmega\allGammaOmega_collected\forHeritability\Yeo_7network_std\RemovingFCeffects\Uni_Yeo7Net_removingFCcov';
path.data='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info';
path.output='F:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info\HOA_20221006';
% % Stats={'AAL','BA','BNA'};
% 
% % for i=1:length(Stats)
%     
%     AE_dataDir=fullfile(path.data,'twins_data_coass_stdYeo7u_forACE_univ_UnivAE.csv');
%     ICC_dataDir=fullfile(path.data,'twins_icc.csv');
%     AIC_pertOutputDir=fullfile(path.data,['ResultsArr.mat']);
%     OutputDir=fullfile(path.data,['ResultsArr_sigTest.csv']);
% ArrageACEresults(AE_dataDir,ICC_dataDir,AIC_pertOutputDir,OutputDir);
% % end

% AE_data=readtable(fullfile(path.data,'twins_data_coass_stdYeo7u_forACE_univ_UnivAE.csv'));
% ACE_data=readtable(fullfile(path.data,'twins_data_coass_stdYeo7u_forACE_univ_UnivACE.csv'));
% CE_data=readtable(fullfile(path.data,'twins_data_coass_stdYeo7u_forACE_univ_UnivCE.csv'));
AE_data=readtable(fullfile(path.data,'twins_data_CBF_HOA_whole_withSmooth_forACE_UnivAE.csv'));
ACE_data=readtable(fullfile(path.data,'twins_data_CBF_HOA_whole_withSmooth_forACE_UnivACE.csv'));
CE_data=readtable(fullfile(path.data,'twins_data_CBF_HOA_whole_withSmooth_forACE_UnivCE.csv'));
ICC_data=readtable(fullfile(path.data,'icc_twin_HOA_whole_withSmooth_forACE.csv'));

ICC_MZ_DZ_sig=ICC_data.p_right<0.05;
ICC_MZ_DZ_sig=double(ICC_MZ_DZ_sig);

p_ACE_AE=AE_data.p_ACE_AE;
p_ACE_CE=AE_data.p_ACE_CE;
p_ACE_E=AE_data.p_ACE_E;

ACE_AIC=AE_data.ACE_AIC;
AE_AIC=AE_data.AE_AIC;
CE_AIC=AE_data.CE_AIC;
E_AIC=AE_data.E_AIC;

AIC_mat=[AE_data.ACE_AIC,AE_data.AE_AIC,AE_data.CE_AIC,AE_data.E_AIC];

[~,Ind_min]=min(AIC_mat,[],2);

% convert data format
if iscell(ACE_data.albound)
ACE_data.albound=cellfun(@str2num,ACE_data.albound,'UniformOutput',0);
 temp_ind=cellfun(@isempty,ACE_data.albound);
 ACE_data.albound(temp_ind==1)=num2cell(0);
 ACE_data.albound=cell2mat(ACE_data.albound);
end
 
if iscell(ACE_data.aubound)
ACE_data.aubound=cellfun(@str2num,ACE_data.aubound,'UniformOutput',0);
temp_ind=cellfun(@isempty,ACE_data.aubound);
 ACE_data.aubound(temp_ind==1)=num2cell(1);
 ACE_data.aubound=cell2mat(ACE_data.aubound);
end

if iscell(ACE_data.clbound)
ACE_data.clbound=cellfun(@str2num,ACE_data.clbound,'UniformOutput',0);
temp_ind=cellfun(@isempty,ACE_data.clbound);
 ACE_data.clbound(temp_ind==1)=num2cell(0);
 ACE_data.clbound=cell2mat(ACE_data.clbound);
end

if iscell(ACE_data.cubound)
ACE_data.cubound=cellfun(@str2num,ACE_data.cubound,'UniformOutput',0);
temp_ind=cellfun(@isempty,ACE_data.cubound);
 ACE_data.cubound(temp_ind==1)=num2cell(1);
 ACE_data.cubound=cell2mat(ACE_data.cubound);
end

if iscell(ACE_data.elbound)
ACE_data.elbound=cellfun(@str2num,ACE_data.elbound,'UniformOutput',0);
temp_ind=cellfun(@isempty,ACE_data.elbound);
 ACE_data.elbound(temp_ind==1)=num2cell(0);
 ACE_data.elbound=cell2mat(ACE_data.elbound);
end

if iscell(ACE_data.eubound)
ACE_data.eubound=cellfun(@str2num,ACE_data.eubound,'UniformOutput',0);
temp_ind=cellfun(@isempty,ACE_data.eubound);
 ACE_data.eubound(temp_ind==1)=num2cell(1);
 ACE_data.eubound=cell2mat(ACE_data.eubound);
end

 %AE
 if iscell(AE_data.albound)
 AE_data.albound=cellfun(@str2num,AE_data.albound,'UniformOutput',0);
 temp_ind=cellfun(@isempty,AE_data.albound);
 AE_data.albound(temp_ind==1)=num2cell(0);
 AE_data.albound=cell2mat(AE_data.albound);
 end
 
 if iscell(AE_data.aubound)
AE_data.aubound=cellfun(@str2num,AE_data.aubound,'UniformOutput',0);
temp_ind=cellfun(@isempty,AE_data.aubound);
 AE_data.aubound(temp_ind==1)=num2cell(1);
 AE_data.aubound=cell2mat(AE_data.aubound);
 end

if iscell(AE_data.elbound)
AE_data.elbound=cellfun(@str2num,AE_data.elbound,'UniformOutput',0);
temp_ind=cellfun(@isempty,AE_data.elbound);
 AE_data.elbound(temp_ind==1)=num2cell(0);
 AE_data.elbound=cell2mat(AE_data.elbound);
end

if iscell(AE_data.eubound)
AE_data.eubound=cellfun(@str2num,AE_data.eubound,'UniformOutput',0);
temp_ind=cellfun(@isempty,AE_data.eubound);
 AE_data.eubound(temp_ind==1)=num2cell(1);
 AE_data.eubound=cell2mat(AE_data.eubound);
end
 % CE
 
 if iscell(CE_data.clbound)
CE_data.clbound=cellfun(@str2num,CE_data.clbound,'UniformOutput',0);
temp_ind=cellfun(@isempty,CE_data.clbound);
 CE_data.clbound(temp_ind==1)=num2cell(0);
 CE_data.clbound=cell2mat(CE_data.clbound);
 end
 
 if iscell(CE_data.cubound)
CE_data.cubound=cellfun(@str2num,CE_data.cubound,'UniformOutput',0);
temp_ind=cellfun(@isempty,CE_data.cubound);
 CE_data.cubound(temp_ind==1)=num2cell(1);
 CE_data.cubound=cell2mat(CE_data.cubound);
 end
 
 if iscell(CE_data.elbound)
CE_data.elbound=cellfun(@str2num,CE_data.elbound,'UniformOutput',0);
temp_ind=cellfun(@isempty,CE_data.elbound);
 CE_data.elbound(temp_ind==1)=num2cell(0);
 CE_data.elbound=cell2mat(CE_data.elbound);
 end
 
 if iscell(CE_data.eubound)
CE_data.eubound=cellfun(@str2num,CE_data.eubound,'UniformOutput',0);
temp_ind=cellfun(@isempty,CE_data.eubound);
 CE_data.eubound(temp_ind==1)=num2cell(1);
 CE_data.eubound=cell2mat(CE_data.eubound);
 end
 
best_model=cell(length(ACE_AIC),1);
best_model_para=zeros(length(ACE_AIC),14);  % 10: MZ_icc 11: DZ_icc 12:p_right 13: PPM
best_model_para(:,10)=ICC_data.r_zyg1;
best_model_para(:,11)=ICC_data.r_zyg2;
best_model_para(:,12)=ICC_data.p_right;
for i=1:length(ACE_AIC)
    
    if ICC_data.p_right(i)<0.05
if p_ACE_AE(i)<0.05 && p_ACE_CE(i)<0.05 && p_ACE_E(i)<0.05
    best_model{i}='ACE';
    
     best_model_para(i,1)=ACE_data.a_2(i);
    
    
   
        best_model_para(i,2)=ACE_data.albound(i);
        best_model_para(i,3)=ACE_data.aubound(i);
        best_model_para(i,4)=ACE_data.c_2(i); 
        best_model_para(i,5)=ACE_data.clbound(i);
      
        best_model_para(i,6)=ACE_data.cubound(i);
       
        best_model_para(i,7)=ACE_data.e_2(i);
        
        best_model_para(i,8)=ACE_data.elbound(i);
         
       best_model_para(i,9)=ACE_data.eubound(i);
        
    
elseif p_ACE_AE(i)>0.05 && p_ACE_CE(i)<0.05 && p_ACE_E(i)<0.05
    best_model{i}='AE';
      
     best_model_para(i,1)=AE_data.a_2(i);
   
        best_model_para(i,2)=AE_data.albound(i);
    
        best_model_para(i,3)=AE_data.aubound(i);
     
        best_model_para(i,7)=AE_data.e_2(i);
        
        best_model_para(i,8)=ACE_data.elbound(i);
        
       best_model_para(i,9)=ACE_data.eubound(i);
        
    
elseif p_ACE_AE(i)<0.05 && p_ACE_CE(i)>0.05 && p_ACE_E(i)<0.05
    best_model{i}='CE';

    
     
        best_model_para(i,4)=CE_data.c_2(i);
      
       
        best_model_para(i,5)=CE_data.clbound(i);
      
        
        best_model_para(i,6)=CE_data.cubound(i);
        
    
        best_model_para(i,7)=CE_data.e_2(i);
        
         
        best_model_para(i,8)=CE_data.elbound(i);
          
        
       best_model_para(i,9)=CE_data.eubound(i);
       
else
    if Ind_min(i)==1
        best_model{i}='ACE';
        
     best_model_para(i,1)=ACE_data.a_2(i);
   
    
        best_model_para(i,2)=ACE_data.albound(i);
    
    
    
        best_model_para(i,3)=ACE_data.aubound(i);
    
      
        best_model_para(i,4)=ACE_data.c_2(i);
      
      
        best_model_para(i,5)=ACE_data.clbound(i);
      
        best_model_para(i,6)=ACE_data.cubound(i);
       
    
         
        best_model_para(i,7)=ACE_data.e_2(i);
        
      
          
        best_model_para(i,8)=ACE_data.elbound(i);
         
       best_model_para(i,9)=ACE_data.eubound(i);
        
    elseif Ind_min(i)==2
        best_model{i}='AE';
         
     best_model_para(i,1)=AE_data.a_2(i);
   
        best_model_para(i,2)=AE_data.albound(i);
    
        best_model_para(i,3)=AE_data.aubound(i);
    
        best_model_para(i,7)=AE_data.e_2(i);
        
        best_model_para(i,8)=ACE_data.elbound(i);
          
       best_model_para(i,9)=ACE_data.eubound(i);
       
    elseif Ind_min(i)==3
        best_model{i}='CE';
       
        best_model_para(i,4)=CE_data.c_2(i);
     
        best_model_para(i,5)=CE_data.clbound(i);
       
        best_model_para(i,6)=CE_data.cubound(i);
        
        best_model_para(i,7)=CE_data.e_2(i);
        
        best_model_para(i,8)=CE_data.elbound(i);
          
       best_model_para(i,9)=CE_data.eubound(i);
        
    elseif Ind_min(i)==4
        best_model{i}='E';
         best_model_para(i,7)=1;
    best_model_para(i,8)=1;
    best_model_para(i,9)=1;
    end
end

    else
      best_model{i}='none';   
    end
% calculate the PPM


if ICC_MZ_DZ_sig(i)==1
    % genetic factors
if strcmp(best_model{i},'ACE')==1
    AIC_diff=CE_AIC(i)-ACE_AIC(i);
    temp_PPMaic=1./(1+exp(-AIC_diff));
    best_model_para(i,13)=temp_PPMaic;
elseif strcmp(best_model{i},'AE')==1

     AIC_diff=E_AIC(i)-AE_AIC(i);
    temp_PPMaic=1./(1+exp(-AIC_diff));
    best_model_para(i,13)=temp_PPMaic;
else
    best_model_para(i,13)=0.5;
end
else
     best_model_para(i,13)=0.5;
end % if ICC_MZ_DZ_sig(i)==1

   % shared environmental factors
if strcmp(best_model{i},'ACE')==1
    AIC_diff=AE_AIC(i)-ACE_AIC(i);
    temp_PPMaic=1./(1+exp(-AIC_diff));
    best_model_para(i,14)=temp_PPMaic;
elseif strcmp(best_model{i},'CE')==1

     AIC_diff=E_AIC(i)-CE_AIC(i);
    temp_PPMaic=1./(1+exp(-AIC_diff));
    best_model_para(i,14)=temp_PPMaic;
else
    best_model_para(i,14)=0.5;
end

    
end



output.modelNames=best_model;
output.best_model_para=best_model_para;
save(fullfile(path.output,'BestModelInfo_3.mat'),'output');