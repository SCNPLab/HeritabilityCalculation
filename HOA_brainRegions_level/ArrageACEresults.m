function ArrageACEresults(AE_dataDir,ICC_dataDir,AIC_pertOutputDir,OutputDir)
%% arrange the ACE results and perform PPM test

% path.data='H:\IPCAS_TWIN\CBF';
% path.output='H:\IPCAS_TWIN\CBF';
% AE_data=readtable(fullfile(dataDir,'twins_data_CBF_AAL_forACE_UnivAE.csv'));
% ICC_data=readtable(fullfile(path.data,'icc_CBF_AAL.csv'));
AE_data=readtable(AE_dataDir);
ICC_data=readtable(ICC_dataDir);
%% determine the best model
AIC_mat=[AE_data.ACE_AIC,AE_data.AE_AIC,AE_data.CE_AIC,AE_data.E_AIC];

[~,Ind_min]=min(AIC_mat,[],2);
AIC_pert=[numel(Ind_min(Ind_min==1))/numel(Ind_min),...
    numel(Ind_min(Ind_min==2))/numel(Ind_min),...
    numel(Ind_min(Ind_min==3))/numel(Ind_min),...
    numel(Ind_min(Ind_min==4))/numel(Ind_min)];

%% calculate the PPM
ICC_MZ_DZ_sig=ICC_data.p_right<0.05;
ICC_MZ_DZ_sig=double(ICC_MZ_DZ_sig);

p_AE_E_sig_con=AE_data.p_AE_E<0.05;
p_AE_E_sig_con=double(p_AE_E_sig_con);
p_AE_E_sig=AE_data.p_AE_E;
p_AE_E_sig(p_AE_E_sig_con==0)=-1;

p_AE_E_sig_ICC_MZ_DZ_sig=ICC_MZ_DZ_sig.*p_AE_E_sig;

minuslogP=zeros(length(p_AE_E_sig_ICC_MZ_DZ_sig),1);

for i=1:length(p_AE_E_sig_ICC_MZ_DZ_sig)
    if p_AE_E_sig_ICC_MZ_DZ_sig(i)>0
        
        minuslogP(i)=-log10(p_AE_E_sig_ICC_MZ_DZ_sig(i));
    else
        minuslogP(i)=0;
    end
end

% significant a_2 after P_AE_E
a_2_sig_AE_E=AE_data.a_2;
p_AE_E_sig_ind=p_AE_E_sig>0;
p_AE_E_sig_ind=double(p_AE_E_sig_ind);
a_2_sig_AE_E=a_2_sig_AE_E.*p_AE_E_sig_ind;
a_2_sig_AE_E(a_2_sig_AE_E==0)=-1;

% significant e_2 after P_AE_E
e_2_sig_AE_E=AE_data.e_2;
% p_AE_E_sig_ind=p_AE_E_sig>0;
% p_AE_E_sig_ind=double(p_AE_E_sig_ind);
e_2_sig_AE_E=e_2_sig_AE_E.*p_AE_E_sig_ind;
e_2_sig_AE_E(e_2_sig_AE_E==0)=-1;

AIC_EminusAE=AE_data.E_AIC-AE_data.AE_AIC;
minuslogP_ind=minuslogP>0;
minuslogP_ind=double(minuslogP_ind);

AIC_Diff_EminusAE=AIC_EminusAE.*minuslogP_ind;

PPMaic=1./(1+exp(-AIC_Diff_EminusAE));

albound=AE_data.albound;
albound=num2cell(albound);
% albound=cell2mat(albound);

aubound=AE_data.aubound;
aubound=num2cell(aubound);
% aubound=strrep(aubound,'NA','10');
% aubound=cell2mat(aubound);

elbound=AE_data.elbound;
elbound=num2cell(elbound);
% elbound=strrep(elbound,'NA','-1');
% elbound=cell2mat(elbound);

eubound=AE_data.eubound;
eubound=num2cell(eubound);
% eubound=strrep(eubound,'NA','-1');
% eubound=cell2mat(eubound);

MZ_icc=ICC_data.r_zyg1;
MZ_icc=num2cell(MZ_icc);

DZ_icc=ICC_data.r_zyg2;
DZ_icc=num2cell(DZ_icc);

p_right=ICC_data.p_right;
p_right=num2cell(p_right);

outputDir=[num2cell(a_2_sig_AE_E),albound,aubound,num2cell(e_2_sig_AE_E),elbound,eubound,MZ_icc,...
          DZ_icc,p_right,num2cell(PPMaic)];
      outputDir=cell2table(outputDir);
      ColNames={'a_2_sig','albound','aubound','e_2_sig','elbound','eubound','MZ_icc','DZ_icc','p_right','PPMaic'};
      RowNames=table2cell(ICC_data(:,1));
      outputDir.Properties.VariableNames=ColNames;
      outputDir.Properties.RowNames=RowNames;
      
%       save(fullfile(path.output,'ResultsArr.mat'),'AIC_pert');
%       writetable(outputDir,fullfile(path.output,'ResultsArr_sigTest.csv'),'delimiter',',','WriteRowNames',1,'WriteVariableNames',1);
%       AIC_pertOutputDir,OutputDir
 save(AIC_pertOutputDir,'AIC_pert');
  writetable(outputDir,OutputDir,'delimiter',',','WriteRowNames',1,'WriteVariableNames',1);
     
