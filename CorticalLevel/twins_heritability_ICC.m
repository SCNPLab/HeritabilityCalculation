 function twins_heritability_ICC(inputfile,outicc,outcorr)
% written in 20161014

%% input
% inputfile = 'J:\TWINIPCAS2015_REST_2mm_scrub_compcor\her_UN\20161014_twins_heritability_new\twins_data_csv_EpA_001.csv';
% outicc = 'J:\TWINIPCAS2015_REST_2mm_scrub_compcor\her_UN\20161014_twins_heritability_new\twins_icc.csv';
% outcorr = 'J:\TWINIPCAS2015_REST_2mm_scrub_compcor\her_UN\20161014_twins_heritability_new\twins_corr.csv';
% inputfile='D:\Owen\CBF\SmoothedData\twins_data_CBF_BNA_forACE.csv';
% outicc='D:\Owen\CBF\SmoothedData\icc_CBF_BNA.csv';
% outcorr='D:\Owen\CBF\SmoothedData\icc_corr_CBF_BNA.csv';
%% read data
fprintf('\tLoading data...\n');
raw_table = readtable(inputfile);



% zyg_mat = cellfun(@str2num, raw_data(:, end));
% beh_img_mat = cellfun(@str2num, raw_data(:, 2:end-1));
zyg_mat=raw_table.zyg;
beh_img_mat=raw_table(:,2:end-1);
beh_img_mat=table2array(beh_img_mat);

num_beh_img = size(beh_img_mat, 2) / 2;
beh_img_mat_1 = beh_img_mat(:, 1:num_beh_img);
beh_img_mat_2 = beh_img_mat(:, (num_beh_img+1):end);

%% Mask outliers
fprintf('\tMask outliers...\n');
beh_img_mat_nan = isnan(beh_img_mat_1 + beh_img_mat_2);
zyg1_good = bsxfun(@and, ~beh_img_mat_nan, zyg_mat == 1);
zyg2_good = bsxfun(@and, ~beh_img_mat_nan, zyg_mat == 2);

%% Calculate ICC
fprintf('\tCalculating Intra-Class Coefficient...\n');

[r, LB, UB, F, df1, df2, p] = arrayfun(@(x) ICC([beh_img_mat_1(zyg1_good(:, x), x), beh_img_mat_2(zyg1_good(:, x), x)], 'C-k'), (1:num_beh_img)');
ICC_mat = [r, LB, UB, F, df1, df2, p];
data1 = r; num1 = sum(zyg1_good)';

[r, LB, UB, F, df1, df2, p] = arrayfun(@(x) ICC([beh_img_mat_1(zyg2_good(:, x), x), beh_img_mat_2(zyg2_good(:, x), x)], 'C-k'), (1:num_beh_img)');
ICC_mat = [ICC_mat, r, LB, UB, F, df1, df2, p];
data2 = r; num2 = sum(zyg2_good)';

% compare r value http://vassarstats.net/rdiff.html
[zval, p_right, p_left, p_both] = fisher_z_compare(data1, data2, num1, num2);

ICC_rst_title = {'r', 'LB', 'UB', 'F', 'df1', 'df2', 'p'};
% ICC_data_strs = [{'behav'}, arrayfun(@(x) sprintf('img_%d', x), 1:(num_beh_img-1), 'UniformOutput', false)]';
ICC_data_strs =arrayfun(@(x) sprintf('img_%d', x), 1:(num_beh_img), 'UniformOutput', false);
ICC_data_strs =ICC_data_strs';
out_title = [{'data', 'zval', 'p_right', 'p_left', 'p_both'}, strcat(ICC_rst_title, '_zyg1'), strcat(ICC_rst_title, '_zyg2')];
out_data = [ICC_data_strs, num2cell([zval, p_right, p_left, p_both, ICC_mat])];
% brant_write_csv(outicc, [out_title; out_data]);
out_data=cell2table(out_data);
out_data.Properties.VariableNames=out_title;
writetable(out_data,outicc);

%% Calculate Pearson's correlation
fprintf('\tCalculating Pearson''s correlation...\n');

[r, p] = arrayfun(@(x) corr(beh_img_mat_1(zyg1_good(:, x), x), beh_img_mat_2(zyg1_good(:, x), x)), (1:num_beh_img)');
corr_mat = [r, p];
data1 = r; num1 = sum(zyg1_good)';

[r, p] = arrayfun(@(x) corr(beh_img_mat_1(zyg2_good(:, x), x), beh_img_mat_2(zyg2_good(:, x), x)), (1:num_beh_img)');
corr_mat = [corr_mat, r, p];
data2 = r; num2 = sum(zyg2_good)';

% compare r value http://vassarstats.net/rdiff.html
[zval, p_right, p_left, p_both] = fisher_z_compare(data1, data2, num1, num2);

corr_title = {'r', 'p'};
% corr_data_strs = [{'behav'}, arrayfun(@(x) sprintf('img_%d', x), 1:(num_beh_img-1), 'UniformOutput', false)]';
corr_data_strs = arrayfun(@(x) sprintf('img_%d', x), 1:(num_beh_img), 'UniformOutput', false)';
out_title = [{'data', 'zval', 'p_right', 'p_left', 'p_both'}, strcat(corr_title, '_zyg1'), strcat(corr_title, '_zyg2')];
out_data = [corr_data_strs, num2cell([zval, p_right, p_left, p_both, corr_mat])];
% brant_write_csv(outcorr, [out_title; out_data]);
out_data=cell2table(out_data);
out_data.Properties.VariableNames=out_title;
writetable(out_data,outcorr);
fprintf('\tFinished!\n');

function [zval, p_right, p_left, p_both] = fisher_z_compare(data1, data2, num1, num2)
% one to one statistics
% http://vassarstats.net/rdiff.html

if any(data1 .^ 2 == 1), data1(data1 .^ 2 == 1) = data1(data1 .^ 2 == 1) * 0.999; end
if any(data2 .^ 2 == 1), data2(data2 .^ 2 == 1) = data2(data2 .^ 2 == 1) * 0.999; end

if any(data1 .^ 2 > 1), data1(data1 .^ 2 > 1) = NaN; end
if any(data2 .^ 2 > 1), data2(data2 .^ 2 > 1) = NaN; end

data1_r2z = 0.5 .* log((1 + data1) ./ (1 - data1));
data2_r2z = 0.5 .* log((1 + data2) ./ (1 - data2));
zval = (data1_r2z - data2_r2z) ./ (sqrt(1 ./ (num1 - 3) + 1 ./ (num2 - 3)));

p_right = 1 - normcdf(zval, 0, 1);
p_left = normcdf(zval, 0, 1);
p_both = 2 * normcdf(-1 * abs(zval), 0, 1);