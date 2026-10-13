function [results, winners] = roi_optimization( ...
    tbv_signal, afni_signal, roi_dir, original_roi_file)

% ROI_OPTIMIZATION
%
% Compare RTP and AFNI signals using:
%   1. Original ROI
%   2. Static ROI optimization
%   3. Dynamic ROI optimization
%
% Inputs
% ------
% tbv_signal          260 x 1 RTP signal
% afni_signal         X x Y x Z x 260 preprocessed AFNI data
% roi_dir             folder containing the 810 candidate ROI masks
% original_roi_file   filename of the original ROI mask
%
% Example:
%
% [results, winners] = roi_optimization( ...
%     tbv_signal, ...
%     afni_signal, ...
%     'amygballs', ...
%     'amygball_21_5_-16+tlrc.BRIK');


n_trs = 260;
win_len = 20;


%% Load candidate ROIs

roi_files = dir(fullfile(roi_dir, '*.BRIK'));

if numel(roi_files) ~= 810
    error('Expected 810 candidate ROIs.')
end


%% Extract AFNI time series from every candidate ROI

candidate_ts = nan(n_trs, numel(roi_files));

for r = 1:numel(roi_files)

    roi_mask = BrikLoad( ...
        fullfile(roi_files(r).folder, roi_files(r).name)) > 0;

    candidate_ts(:,r) = roi_time_series(afni_signal, roi_mask);

end


%% Original ROI

original_mask = BrikLoad( ...
    fullfile(roi_dir, original_roi_file)) > 0;

original_ts = roi_time_series(afni_signal, original_mask);


%% Define task blocks
%
% 1 = Rest
% 2 = Happy
% 3 = Count

block_order = repmat([1 2 3], 1, ceil(260/20/3));
block_order = block_order(1:13);

happy_blocks = find(block_order == 2);
rest_blocks  = find(block_order == 1);

happy_mask = false(n_trs,1);
rest_mask  = false(n_trs,1);

for b = 1:13

    idx = (b-1)*20 + 1 : b*20;

    if block_order(b) == 2
        happy_mask(idx) = true;

    elseif block_order(b) == 1
        rest_mask(idx) = true;
    end

end


%% Original correlations

corr_original_full = corr( ...
    tbv_signal, original_ts, ...
    'rows', 'complete');

corr_original_happy = corr( ...
    tbv_signal(happy_mask), ...
    original_ts(happy_mask), ...
    'rows', 'complete');

corr_original_rest = corr( ...
    tbv_signal(rest_mask), ...
    original_ts(rest_mask), ...
    'rows', 'complete');


%% Static optimization
%
% Select the candidate ROI with the highest RTP-AFNI correlation.

[corr_static_full, static_full_idx] = ...
    best_roi(tbv_signal, candidate_ts);

[corr_static_happy, static_happy_idx] = ...
    best_roi( ...
        tbv_signal(happy_mask), ...
        candidate_ts(happy_mask,:));

[corr_static_rest, static_rest_idx] = ...
    best_roi( ...
        tbv_signal(rest_mask), ...
        candidate_ts(rest_mask,:));


%% Dynamic optimization: full run
%
% Twenty-TR windows move forward one TR at a time.
% This gives 241 overlapping windows in a 260-TR run.
%
% The winning AFNI time series from each window is added to a cumulative
% signal. Overlapping values are averaged before the final correlation.

n_windows = n_trs - win_len + 1;

dynamic_sum = zeros(n_trs,1);
weights = zeros(n_trs,1);

dynamic_full_idx = nan(n_windows,1);

for w = 1:n_windows

    idx = w : w + win_len - 1;

    [~, best_idx] = best_roi( ...
        tbv_signal(idx), ...
        candidate_ts(idx,:));

    dynamic_sum(idx) = ...
        dynamic_sum(idx) + candidate_ts(idx,best_idx);

    weights(idx) = weights(idx) + 1;

    dynamic_full_idx(w) = best_idx;

end

dynamic_full_ts = dynamic_sum ./ weights;

corr_dynamic_full = corr( ...
    tbv_signal, dynamic_full_ts, ...
    'rows', 'complete');


%% Dynamic optimization: Happy
%
% Select the best ROI separately within each actual 20-TR Happy block.

[corr_dynamic_happy, dynamic_happy_idx] = ...
    dynamic_blocks( ...
        tbv_signal, ...
        candidate_ts, ...
        happy_blocks, ...
        win_len);


%% Dynamic optimization: Rest

[corr_dynamic_rest, dynamic_rest_idx] = ...
    dynamic_blocks( ...
        tbv_signal, ...
        candidate_ts, ...
        rest_blocks, ...
        win_len);


%% Results

Condition = ["FullRun"; "Happy"; "Rest"];

Original = [
    corr_original_full
    corr_original_happy
    corr_original_rest
];

Static = [
    corr_static_full
    corr_static_happy
    corr_static_rest
];

Dynamic = [
    corr_dynamic_full
    corr_dynamic_happy
    corr_dynamic_rest
];

results = table( ...
    Condition, ...
    Original, ...
    Static, ...
    Dynamic);


%% Save winning ROI filenames

winners.Static.FullRun = ...
    roi_files(static_full_idx).name;

winners.Static.Happy = ...
    roi_files(static_happy_idx).name;

winners.Static.Rest = ...
    roi_files(static_rest_idx).name;


winners.Dynamic.FullRun = ...
    {roi_files(dynamic_full_idx).name}';

winners.Dynamic.Happy = ...
    {roi_files(dynamic_happy_idx).name}';

winners.Dynamic.Rest = ...
    {roi_files(dynamic_rest_idx).name}';


end


%% ========================================================================
% Mean AFNI time series within one ROI
% ========================================================================

function ts = roi_time_series(data, mask)

n = size(data,4);

voxels = data(repmat(mask, [1 1 1 n]));

ts = mean( ...
    reshape(voxels, [], n), ...
    1)';

end


%% ========================================================================
% Find candidate ROI with highest correlation
% ========================================================================

function [best_r, best_idx] = best_roi(rtp, candidate_ts)

best_r = -Inf;
best_idx = NaN;

for r = 1:size(candidate_ts,2)

    this_r = corr( ...
        rtp, ...
        candidate_ts(:,r), ...
        'rows', 'complete');

    if this_r > best_r

        best_r = this_r;
        best_idx = r;

    end

end

end


%% ========================================================================
% Dynamic optimization within actual task blocks
% ========================================================================

function [final_r, winning_idx] = ...
    dynamic_blocks(rtp, candidate_ts, blocks, block_length)

rtp_concat = [];
afni_concat = [];

winning_idx = nan(numel(blocks),1);

for b = 1:numel(blocks)

    block = blocks(b);

    idx = ...
        (block-1)*block_length + 1 : ...
        block*block_length;

    [~, best_idx] = best_roi( ...
        rtp(idx), ...
        candidate_ts(idx,:));

    rtp_concat = [
        rtp_concat
        rtp(idx)
    ];

    afni_concat = [
        afni_concat
        candidate_ts(idx,best_idx)
    ];

    winning_idx(b) = best_idx;

end


final_r = corr( ...
    rtp_concat, ...
    afni_concat, ...
    'rows', 'complete');

end
