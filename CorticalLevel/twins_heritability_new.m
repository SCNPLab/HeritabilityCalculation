function twins_heritability_new(inputfile,outrandfn,outrstfn)
% written in 20161014
% path.root='H:\IPCAS_TWIN\CBF\SmoothedData';
%% input
cov_strs = {'gender','age','FD_power'};
zyg_str = 'gyz';
% beh_str = 'ZRE_1'; % can only be one
% data_str = 'EpA_*'; % all matched columns
%inputfile = 'J:\TWINIPCAS2015_REST_2mm_scrub_compcor\her_UN\20161014_twins_heritability_new\202twins_EP_A_selfconnectionHz.csv';
% outrandfn = 'J:\TWINIPCAS2015_REST_2mm_scrub_compcor\her_UN\20161014_twins_heritability_new\randomize_info.csv';
% outrstfn = 'J:\TWINIPCAS2015_REST_2mm_scrub_compcor\her_UN\20161014_twins_heritability_new\twins_data_csv_EpA_001.csv';
% inputfile=fullfile(path.root,'CBF_BA_forACE.csv');  %%  AAL BA BNA
% outrandfn=fullfile(path.root,'randomize_info_CBF_BA.csv');
% outrstfn=fullfile(path.root,'twins_data_CBF_BA_forACE.csv');

% raw_table = brant_read_csv(inputfile);
raw_table=readtable(inputfile);
raw_title = raw_table.Properties.VariableNames;
% raw_data = raw_table(2:end, :);

%% check for subjects' IDs (missing or unexpected)
fprintf('\tChecking for IDs (missing or unexpected)...\n');
% find columns of ID. twin_ids: unique string of each twin-pair
% subj_list: raw list of subject ids
% ID_ind = find(strcmpi(raw_title, 'ID'));
% subj_list = raw_data(:, ID_ind); %#ok<FNDSB>
subj_list=raw_table.ID;
subj_list=num2str(subj_list);
subj_list=cellstr(subj_list);
twin_ids = unique(cellfun(@(x) x(1:end-2), subj_list, 'UniformOutput', false), 'stable');
subj_list_full = [strcat(twin_ids, '01'); strcat(twin_ids, '02')];

subj_missing = setdiff(subj_list_full, subj_list);
if ~isempty(subj_missing)
    error([sprintf('%s\n', subj_missing{:}),...
           sprintf('Above subjects are missing!\n')]);
end

subj_werid = setdiff(subj_list, subj_list_full);
if ~isempty(subj_werid)
    error([sprintf('%s\n', subj_werid{:}),...
           sprintf('Above subjects are unexpected!\n')]);
end

%% get zyg data
fprintf('\tGeting zyg information...\n');
% zyg_ind = find(strcmpi(zyg_str, raw_title));
% assert(numel(zyg_ind) == 1);
% zyg_mat = cellfun(@str2num, raw_data(:, zyg_ind));
zyg_mat=raw_table.(zyg_str);
%% get covariates
fprintf('\tGeting covariates...\n');
% cov_inds = cellfun(@(x) find(strcmpi(x, raw_title)), cov_strs);
% cov_mat = cellfun(@str2num, raw_data(:, cov_inds));
     cov_mat=[];
    for cov_strs_ord=1:length(cov_strs)
       cov_mat = [cov_mat,raw_table.(cov_strs{cov_strs_ord})];
    end
cov_mat_ones = [ones(size(cov_mat, 1), 1), cov_mat];

%% get behavior data
% fprintf('\tGeting behavior data...\n');
% beh_ind = find(strcmpi(beh_str, raw_title));
% assert(numel(beh_ind) == 1);
% beh_mat = cellfun(@str2num, raw_data(:, beh_ind));

%% get image data
fprintf('\tGeting image data...\n');
% data_match_ind_tmp = cellfun(@(x) regexp(x, data_str), raw_title, 'UniformOutput', false);
% data_match_ind = ~cellfun(@isempty, data_match_ind_tmp);
% img_mat = cellfun(@str2num, raw_data(:, data_match_ind));
img_mat=raw_table(:,6:end);
img_mat=table2array(img_mat);

%% 1. regress out covariates
% regress with all subject at one time.
fprintf('\tRegressing out covariates for behavior and image data...\n');
% beh_mat_res = beh_mat - cov_mat_ones * (cov_mat_ones \ beh_mat);
img_mat_res = img_mat - cov_mat_ones * (cov_mat_ones \ img_mat);

%% 2. normalise for image data and behavior data (rawdata - mean)/std
fprintf('\tNormalize for behavior and image data... (using zscore, (rawdata - mean)/std)\n');
% beh_mat_res_nor = zscore(beh_mat_res); %(beh_mat_res - mean(beh_mat_res)) ./ std(beh_mat_res);
img_mat_res_nor = zscore(img_mat_res); %bsxfun(@rdivide, bsxfun(@minus, img_mat_res, mean(img_mat_res)), std(img_mat_res));

%% 3. randomize for each pair of twins
fprintf('\tDoing randomization for each pair of twins...\n');
rng('default');
rand_twin = zeros(numel(twin_ids), 2);
for m = 1:numel(twin_ids)
    rand_twin(m, :) = randperm(2);
end
rand_twin_strs = arrayfun(@(x, y) sprintf('%s%02d', x{1}, y), [twin_ids, twin_ids], rand_twin, 'UniformOutput', false);
rand_twin_inds = cellfun(@(x) find(strcmpi(subj_list, x)), rand_twin_strs);

% beh_mat_res_nor_rand_1 = beh_mat_res_nor(rand_twin_inds(:, 1), :);
% beh_mat_res_nor_rand_2 = beh_mat_res_nor(rand_twin_inds(:, 2), :);
img_mat_res_nor_rand_1 = img_mat_res_nor(rand_twin_inds(:, 1), :);
img_mat_res_nor_rand_2 = img_mat_res_nor(rand_twin_inds(:, 2), :);

% new order of zyg
zyg_mat_rand_1 = zyg_mat(rand_twin_inds(:, 1), :);
zyg_mat_rand_2 = zyg_mat(rand_twin_inds(:, 2), :);
assert(isequal(zyg_mat_rand_1, zyg_mat_rand_2));
zyg_mat_rand = zyg_mat_rand_1;

% %% make ourliers NaN seperately
% beh_mat_res_nor_rand_1 = remove_outliers_sep(beh_mat_res_nor_rand_1);
% beh_mat_res_nor_rand_2 = remove_outliers_sep(beh_mat_res_nor_rand_2);
% img_mat_res_nor_rand_1 = remove_outliers_sep(img_mat_res_nor_rand_1);
% img_mat_res_nor_rand_2 = remove_outliers_sep(img_mat_res_nor_rand_2);

%% make ourliers NaN as whole
fprintf('\tMake ourliers NaN as whole...\n');
% [beh_mat_res_nor_rand_1, beh_mat_res_nor_rand_2] = remove_outliers_whole([beh_mat_res_nor_rand_1, beh_mat_res_nor_rand_2]);
[img_mat_res_nor_rand_1, img_mat_res_nor_rand_2] = remove_outliers_whole([img_mat_res_nor_rand_1, img_mat_res_nor_rand_2]);

%% output randomization information to csv
fprintf('\tOutput randomization information to %s...\n', outrandfn);
out_rand_title = {'ID', 'twin1', 'twin2', 'zyg'};
out_rand_data = [twin_ids, rand_twin_strs, num2cell(zyg_mat_rand)];
% brant_write_csv(outrandfn, [out_rand_title; out_rand_data]);
out_rand_data=cell2table(out_rand_data);
out_rand_data.Properties.VariableNames=out_rand_title;
writetable(out_rand_data,outrandfn);

%% output regressed, standardized and randomized and outlier removed to csv
fprintf('\tOutput regressed, standardized, randomized and outlier removed data to %s...\n', outrstfn);
num_img_data = size(img_mat_res_nor, 2);
twin1_img_title = arrayfun(@(x) sprintf('twin1_%d', x), 1:num_img_data, 'UniformOutput', false);
twin2_img_title = arrayfun(@(x) sprintf('twin2_%d', x), 1:num_img_data, 'UniformOutput', false);

%out_title = [{'ID'}, {'twin1_behav'}, twin1_img_title, {'twin2_behav'}, twin2_img_title, {'zyg'}];
%out_data = [twin_ids, num2cell([beh_mat_res_nor_rand_1, img_mat_res_nor_rand_1, beh_mat_res_nor_rand_2, img_mat_res_nor_rand_2, zyg_mat_rand])];
out_title = [{'ID'}, twin1_img_title,  twin2_img_title, {'zyg'}];
out_data = [twin_ids,num2cell([img_mat_res_nor_rand_1, img_mat_res_nor_rand_2, zyg_mat_rand])];

% brant_write_csv(outrstfn, [out_title; out_data]);
out_data=cell2table(out_data);
out_data.Properties.VariableNames=out_title;
writetable(out_data,outrstfn);
fprintf('\tFinished!\n');

function data_mat = remove_outliers_sep(data_mat)
%the function of exlcding the outliers based on all of the subjects: 3 IQR
for m = 1:size(data_mat, 2)
    tmp = data_mat(:, m);
    med_mat = median(tmp);
    Q1 = median(tmp(tmp < med_mat));
    Q3 = median(tmp(tmp > med_mat));
    iqr_mat = iqr(tmp);

    %criterion of selecting bad voxel
    bad_ind = (tmp > (Q3 + 3 * iqr_mat)) | (tmp < (Q1 - 3 * iqr_mat));
    %replace the value of bad voxel with "NaN"
    data_mat(bad_ind, m) = NaN;
end

function [data_mat_out1, data_mat_out2] = remove_outliers_whole(data_mat)
%the function of exlcding the outliers based on all of the subjects: 3 IQR
size_in = size(data_mat);
data_mat_tmp = [data_mat(:, 1:size_in(2)/2); data_mat(:, size_in(2)/2+1:end)];

for m = 1:size(data_mat_tmp, 2)
    tmp = data_mat_tmp(:, m);
    med_mat = median(tmp);
    Q1 = median(tmp(tmp < med_mat));
    Q3 = median(tmp(tmp > med_mat));
    iqr_mat = iqr(tmp);

    %criterion of selecting bad voxel
    bad_ind = (tmp > (Q3 + 3 * iqr_mat)) | (tmp < (Q1 - 3 * iqr_mat));
    %replace the value of bad voxel with "NaN"
    data_mat_tmp(bad_ind, m) = NaN;
end

data_mat_out1 = data_mat_tmp(1:size_in(1), :);
data_mat_out2 = data_mat_tmp(size_in(1)+1:end, :);