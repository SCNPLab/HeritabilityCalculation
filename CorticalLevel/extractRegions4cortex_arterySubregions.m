clear all;
path.artery='D:\NeuroImagingSoftwares\BrainAtlas\mni_vascular_territories.nii';
path.HOA='H:\IPCAS_TWIN\CBF\info\Atlases\Reslice_Atlases';
path.label='H:\IPCAS_TWIN\CBF\info\Info';
path.output='H:\IPCAS_TWIN\CBF\SmoothedData\20220620\Info\average4Cortex_HOV\reference_images\DivBilateral_sliced';
if ~exist(path.output,'dir')
    mkdir(path.output)
end


HOA_V=spm_vol(fullfile(path.HOA,'HOC_cortical_subcortical.nii'));
HOA_img=spm_read_vols(HOA_V);

HOV_remove_brainstem=HOA_img;
HOV_remove_brainstem(HOV_remove_brainstem==97)=0;
HOV_remove_brainstem(HOV_remove_brainstem==98)=0;
HOV_remove_brainstem_V=HOA_V;
HOV_remove_brainstem_V.fname=fullfile(path.HOA,'HOC_cortical_subcortical_removeBrainstem.nii');
spm_write_vol(HOV_remove_brainstem_V,HOV_remove_brainstem);


HOA_label=readtable(fullfile(path.label,'HOA_labels_cortexLabels.xlsx'));
HOA_cortex_label=HOA_label.Cortex_label;
unique_label=unique(HOA_cortex_label);
unique_label(unique_label==0)=[];
unique_label_names={'L-Frontal','R-Frontal','L-Insula','R-Insula','L-Temporal','R-Temporal','L-Parietal','R-Parietal','L-Occipital','R-Occipital','L-Cingulate','R-Cingulate','L-Subcortical','R-Subcortical'};

Whole_HOA_cortex=zeros(size(HOA_img));
Whole_HOA_cortex_V=HOA_V;

for i=13:length(unique_label)
    
    index_label=find(HOA_cortex_label==unique_label(i));
    HOA_output=zeros(size(HOA_img));
    
    for k = 1:length(index_label)
       
        HOA_output(HOA_img==index_label(k))=1;
        Whole_HOA_cortex(HOA_img==index_label(k))=i;
    end
    
    HOA_output_V=HOA_V;
    HOA_output_V.fname=fullfile(path.output,['HOA_' unique_label_names{i} '.nii']);
     spm_write_vol(HOA_output_V,HOA_output);
end
Whole_HOA_cortex_V.fname=fullfile(path.output,'HOA_diff_cortices.nii');
%     spm_write_vol(Whole_HOA_cortex_V,Whole_HOA_cortex);

% ACA MCA PCA
artery_V=spm_vol(fullfile(path.artery,'Reslice_mni_vascular_territories.nii'));
artery_img=spm_read_vols(artery_V);
unique_artery=unique(artery_img);
unique_artery(unique_artery==0)=[];

remain_cortex=zeros(size(artery_img));
% % remain_cortex(remain_cortex==4)=0;
% % remain_cortex(remain_cortex==5)=0;
% % remain_cortex(remain_cortex==9)=0;
% % remain_cortex(remain_cortex==10)=0;
remain_cortex_V=artery_V;
remain_cortex_V.fname=fullfile(path.output,'artery_remain_cortex.nii');


for ii=1:length(unique_artery)
    artery_output=zeros(size(artery_img));
    artery_output(artery_img==unique_artery(ii))=1;
    
    artery_output_V=artery_V;
    artery_output_V.fname=fullfile(path.output,['artery_' num2str(ii) '.nii']);
    if ii~=4 && ii~=5 && ii~=9 && ii~=10
      remain_cortex(artery_img==unique_artery(ii))=ii; 
    end
%     spm_write_vol(artery_output_V,artery_output);
end
        
% spm_write_vol(remain_cortex_V,remain_cortex);