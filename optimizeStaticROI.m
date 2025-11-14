function [bestROI, bestCorr, allCorrs] = optimizeStaticROI(tbv_signal, afni_signal, roi_files)
%OPTIMIZESTATICROI Selects the ROI with highest mean correlation to RTP
%
% INPUTS:
%   tbv_signal   - numeric vector [TRs x 1], RTP signal for the run
%   afni_signal  - 4-D array [x y z t], AFNI run data
%   roi_files    - cell array of strings, paths to ROI .BRIK files
%
% OUTPUTS:
%   bestROI      - string, filename of ROI with highest correlation
%   bestCorr     - scalar, correlation of best ROI
%   allCorrs     - table, columns: ROI (file name), Corr (correlation)
%
% Example usage:
%   [bestROI, bestCorr, allCorrs] = optimizeStaticROI(tbv, afni, roi_files);

nROIs = numel(roi_files);
allCorrs = table('Size',[nROIs 2],'VariableTypes',{'string','double'},...
                 'VariableNames',{'ROI','Corr'});

bestCorr = -Inf;
bestROI = '';

for r = 1:nROIs
    roi_mask = BrikLoad(roi_files{r}) > 0;           % load ROI mask
    roi_voxels = afni_signal(repmat(roi_mask,[1 1 1 size(afni_signal,4)]));
    roi_ts = mean(reshape(roi_voxels, [], size(afni_signal,4)), 1)';  % mean across voxels
    r_val = corr(tbv_signal, roi_ts, 'rows', 'complete');
    
    allCorrs.ROI(r) = string(roi_files{r});
    allCorrs.Corr(r) = r_val;
    
    if r_val > bestCorr
        bestCorr = r_val;
        bestROI = roi_files{r};
    end
end

end
