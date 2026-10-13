% 04_procedural_null_simulation.m
%
% Procedural-null simulation for RTP-AFNI ROI optimization.
%
% The simulation asks how much correlation gain can arise from searching
% across candidate ROIs when the RTP and AFNI signals contain shared signal
% plus independent noise.



%% Parameters

T = 260;
n_simulations = 100;

signal_weight = 0.7;

lambda = 0.2;
boost_lambda = 2;
top_n = 500;

fwhm_mm = 4;
voxel_size = 2;

sigma = (fwhm_mm / voxel_size) / 2.355;

win_len = 20;


%% Input files

% This dataset is used only to define the dimensions of the simulated volume.
sample_volume = BrikLoad('002_1_corr_roi1_rtp_BL+tlrc.BRIK');

vol_size = size(sample_volume);

if length(vol_size) > 3
    vol_size = vol_size(1:3);
end


roi_dirs = {
    'amygballs'
    'hipsballs'
};

target_names = {
    'Amygdala'
    'HIPS'
};

original_roi_files = {
    'amygballs/amygball_21_5_-16+tlrc.BRIK'
    'hipsballs/hipsball_42_48_48+tlrc.BRIK'
};


%% Load candidate masks

candidate_masks = cell(2,1);
original_masks = cell(2,1);

for target = 1:2

    roi_files = dir( ...
        fullfile(roi_dirs{target}, '*.BRIK'));

    if numel(roi_files) ~= 810
        error('Expected 810 candidate ROIs.')
    end

    candidate_masks{target} = cell(810,1);

    for r = 1:810

        candidate_masks{target}{r} = ...
            BrikLoad( ...
                fullfile(roi_files(r).folder, roi_files(r).name)) > 0;

    end

    original_masks{target} = ...
        BrikLoad(original_roi_files{target}) > 0;

end


%% ========================================================================
% FULL-RUN SIMULATIONS
%
% 100 independent simulations are generated for each target.
% ========================================================================

fullrun_results = table();


for target = 1:2

    original_mask = original_masks{target};
    masks = candidate_masks{target};

    for simulation = 1:n_simulations

        fprintf( ...
            'Full run: %s simulation %d of %d\n', ...
            target_names{target}, ...
            simulation, ...
            n_simulations);


        %% Shared signal and simulated RTP

        true_signal = randn(1,T);

        rtp_signal = ...
            true_signal * 0.6 + ...
            randn(1,T);


        %% Spatial signal weights

        weights = ...
            poissrnd(lambda, vol_size);

        sorted_weights = ...
            sort(weights(:), 'descend');

        top_value = ...
            sorted_weights(top_n);


        % Give voxels within the target ROI stronger signal weights

        boost_mask = ...
            original_mask & ...
            weights < top_value;

        weights(boost_mask) = ...
            top_value + ...
            poissrnd( ...
                boost_lambda, ...
                sum(boost_mask(:)), ...
                1);


        weights = ...
            weights / max(weights(:));


        %% Generate simulated AFNI data

        afni_sim = ...
            zeros([vol_size T]);


        for t = 1:T

            signal_vol = ...
                imgaussfilt3( ...
                    signal_weight * ...
                    true_signal(t) * ...
                    weights, ...
                    sigma);

            noise_vol = ...
                randn(vol_size);

            smooth_noise = ...
                imgaussfilt3( ...
                    noise_vol, ...
                    sigma);

            afni_sim(:,:,:,t) = ...
                signal_vol + smooth_noise;

        end


        %% Extract original ROI signal

        original_ts = ...
            roi_time_series( ...
                afni_sim, ...
                original_mask);


        corr_original = ...
            corr( ...
                rtp_signal', ...
                original_ts);


        %% Extract all candidate ROI signals

        candidate_ts = ...
            nan(T,810);


        for r = 1:810

            candidate_ts(:,r) = ...
                roi_time_series( ...
                    afni_sim, ...
                    masks{r});

        end


        %% Static optimization

        best_static = -Inf;


        for r = 1:810

            this_r = ...
                corr( ...
                    rtp_signal', ...
                    candidate_ts(:,r));

            if this_r > best_static
                best_static = this_r;
            end

        end


        %% Dynamic optimization
        %
        % 20-TR windows move forward one TR at a time.
        % The winning signal from each window contributes to the
        % reconstructed 260-TR dynamic signal.

        dynamic_sum = zeros(T,1);
        dynamic_weights = zeros(T,1);


        for w = 1:(T - win_len + 1)

            idx = ...
                w : w + win_len - 1;

            best_r = -Inf;
            best_ts = [];


            for r = 1:810

                this_r = ...
                    corr( ...
                        rtp_signal(idx)', ...
                        candidate_ts(idx,r));

                if this_r > best_r

                    best_r = this_r;
                    best_ts = candidate_ts(idx,r);

                end

            end


            dynamic_sum(idx) = ...
                dynamic_sum(idx) + best_ts;

            dynamic_weights(idx) = ...
                dynamic_weights(idx) + 1;

        end


        dynamic_ts = ...
            dynamic_sum ./ dynamic_weights;


        corr_dynamic = ...
            corr( ...
                rtp_signal', ...
                dynamic_ts);


        %% Save result

        new_row = table( ...
            string(target_names{target}), ...
            simulation, ...
            corr_original, ...
            best_static, ...
            corr_dynamic, ...
            'VariableNames', { ...
                'Target', ...
                'Simulation', ...
                'Corr_Original', ...
                'Corr_Static', ...
                'Corr_Dynamic' ...
            });


        fullrun_results = [
            fullrun_results
            new_row
        ];

    end

end


writetable( ...
    fullrun_results, ...
    'procedural_null_fullrun.csv');


%% ========================================================================
% CONDITION-SPECIFIC SIMULATIONS
%
% This is a separate simulation batch from the full-run analysis.
%
% Each simulation produces both Happy and Rest estimates.
% ========================================================================


% Task structure:
% Rest, Happy, Count repeated across 13 20-TR blocks

block_order = ...
    repmat([1 2 3], 1, ceil(13/3));

block_order = ...
    block_order(1:13);


happy_blocks = ...
    find(block_order == 2);

rest_blocks = ...
    find(block_order == 1);


happy_mask = false(T,1);
rest_mask = false(T,1);


for b = 1:13

    idx = ...
        (b-1)*20 + 1 : b*20;

    if block_order(b) == 2

        happy_mask(idx) = true;

    elseif block_order(b) == 1

        rest_mask(idx) = true;

    end

end


condition_results = table();


for target = 1:2

    original_mask = original_masks{target};
    masks = candidate_masks{target};


    for simulation = 1:n_simulations

        fprintf( ...
            'Conditions: %s simulation %d of %d\n', ...
            target_names{target}, ...
            simulation, ...
            n_simulations);


        %% Shared signal and simulated RTP

        true_signal = randn(1,T);

        rtp_signal = ...
            true_signal * 0.6 + ...
            randn(1,T);


        %% Spatial signal weights

        weights = ...
            poissrnd(lambda, vol_size);

        sorted_weights = ...
            sort(weights(:), 'descend');

        top_value = ...
            sorted_weights(top_n);


        boost_mask = ...
            original_mask & ...
            weights < top_value;

        weights(boost_mask) = ...
            top_value + ...
            poissrnd( ...
                boost_lambda, ...
                sum(boost_mask(:)), ...
                1);


        weights = ...
            weights / max(weights(:));


        %% Generate simulated AFNI data

        afni_sim = ...
            zeros([vol_size T]);


        for t = 1:T

            signal_vol = ...
                imgaussfilt3( ...
                    signal_weight * ...
                    true_signal(t) * ...
                    weights, ...
                    sigma);

            noise_vol = ...
                randn(vol_size);

            smooth_noise = ...
                imgaussfilt3( ...
                    noise_vol, ...
                    sigma);

            afni_sim(:,:,:,t) = ...
                signal_vol + smooth_noise;

        end


        %% Extract original and candidate ROI time series

        original_ts = ...
            roi_time_series( ...
                afni_sim, ...
                original_mask);


        candidate_ts = ...
            nan(T,810);


        for r = 1:810

            candidate_ts(:,r) = ...
                roi_time_series( ...
                    afni_sim, ...
                    masks{r});

        end


        %% ---------------------------------------------------------------
        % HAPPY
        % ---------------------------------------------------------------

        corr_original_happy = ...
            corr( ...
                rtp_signal(happy_mask)', ...
                original_ts(happy_mask));


        % Static Happy

        best_static_happy = -Inf;


        for r = 1:810

            this_r = ...
                corr( ...
                    rtp_signal(happy_mask)', ...
                    candidate_ts(happy_mask,r));

            if this_r > best_static_happy
                best_static_happy = this_r;
            end

        end


        % Dynamic Happy:
        % optimize separately within each actual 20-TR Happy block

        rtp_happy = [];
        dynamic_happy = [];


        for b = 1:numel(happy_blocks)

            block = happy_blocks(b);

            idx = ...
                (block-1)*20 + 1 : block*20;

            best_r = -Inf;
            best_ts = [];


            for r = 1:810

                this_r = ...
                    corr( ...
                        rtp_signal(idx)', ...
                        candidate_ts(idx,r));

                if this_r > best_r

                    best_r = this_r;
                    best_ts = candidate_ts(idx,r);

                end

            end


            rtp_happy = [
                rtp_happy
                rtp_signal(idx)'
            ];

            dynamic_happy = [
                dynamic_happy
                best_ts
            ];

        end


        corr_dynamic_happy = ...
            corr( ...
                rtp_happy, ...
                dynamic_happy);


        %% ---------------------------------------------------------------
        % REST
        % ---------------------------------------------------------------

        corr_original_rest = ...
            corr( ...
                rtp_signal(rest_mask)', ...
                original_ts(rest_mask));


        % Static Rest

        best_static_rest = -Inf;


        for r = 1:810

            this_r = ...
                corr( ...
                    rtp_signal(rest_mask)', ...
                    candidate_ts(rest_mask,r));

            if this_r > best_static_rest
                best_static_rest = this_r;
            end

        end


        % Dynamic Rest:
        % optimize separately within each actual 20-TR Rest block

        rtp_rest = [];
        dynamic_rest = [];


        for b = 1:numel(rest_blocks)

            block = rest_blocks(b);

            idx = ...
                (block-1)*20 + 1 : block*20;

            best_r = -Inf;
            best_ts = [];


            for r = 1:810

                this_r = ...
                    corr( ...
                        rtp_signal(idx)', ...
                        candidate_ts(idx,r));

                if this_r > best_r

                    best_r = this_r;
                    best_ts = candidate_ts(idx,r);

                end

            end


            rtp_rest = [
                rtp_rest
                rtp_signal(idx)'
            ];

            dynamic_rest = [
                dynamic_rest
                best_ts
            ];

        end


        corr_dynamic_rest = ...
            corr( ...
                rtp_rest, ...
                dynamic_rest);


        %% Save Happy and Rest results

        happy_row = table( ...
            string(target_names{target}), ...
            simulation, ...
            "Happy", ...
            corr_original_happy, ...
            best_static_happy, ...
            corr_dynamic_happy, ...
            'VariableNames', { ...
                'Target', ...
                'Simulation', ...
                'Condition', ...
                'Corr_Original', ...
                'Corr_Static', ...
                'Corr_Dynamic' ...
            });


        rest_row = table( ...
            string(target_names{target}), ...
            simulation, ...
            "Rest", ...
            corr_original_rest, ...
            best_static_rest, ...
            corr_dynamic_rest, ...
            'VariableNames', { ...
                'Target', ...
                'Simulation', ...
                'Condition', ...
                'Corr_Original', ...
                'Corr_Static', ...
                'Corr_Dynamic' ...
            });


        condition_results = [
            condition_results
            happy_row
            rest_row
        ];

    end

end


writetable( ...
    condition_results, ...
    'procedural_null_conditions.csv');


%% ========================================================================
% Helper function
% ========================================================================

function ts = roi_time_series(data, mask)

    n_trs = size(data,4);

    voxels = ...
        data(repmat(mask, [1 1 1 n_trs]));

    ts = ...
        mean( ...
            reshape(voxels, [], n_trs), ...
            1)';

end
