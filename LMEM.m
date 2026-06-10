%% Network Pair Linear Mixed Effect Model 
%% This script loads subject-level effective connectivity matrices, computes mean connectivity values for predefined between-network masks, and performs group comparisons between healthy controls (HC) and opioid use disorder (OUD) participants.
%%
%% For each network pair, the script:

%% 1. Extracts mean connectivity values from the network mask.
%% 2. Compares groups using a linear mixed-effects model
%%    (with automatic fallback to a standard linear model if required).
%% 3. Computes descriptive statistics, partial eta squared, and
%%    Hedges' g effect sizes.
%% 4. Applies Benjamini-Hochberg false discovery rate (BH-FDR)
%%    correction across all network comparisons.
%% 5. Exports a summary table of statistical results to CSV.
%%
%% Inputs:
%%    - Subject-level effective connectivity matrices (*.csv)
%%    - templateMat: network-pair masks
%%    - TmpNetName: network-pair labels
%%
%% Outputs:
%%    - Network_LMEM_Stats_<Task>.csv
%%
%% Author: Brianna Austin (Imperial College London)
%% Date: 06/05/2026

clear; clc;

%% 1. Firstly, set your paths
feature('numthreads', 1);

baseDir = '/Users/briannaaustin/Downloads/MSc_BriannaA-2';
outputpath = fullfile(baseDir, 'Outputs');

load('templateMatAndNetName.mat');

% Toggle between 'CueData' and 'MIDData'
taskType = 'MIDData';

baseloc_hc = fullfile('/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)', taskType, 'HC_MID');
baseloc_pt = fullfile('/Users/briannaaustin/Desktop/lsngc(2)/EC_Brianna(2)', taskType, 'Patients_MID');

%% 2. Load your HC and PT files from their directories
FilesHC = dir(fullfile(baseloc_hc, '*_Aff.csv'));
FilesPT = dir(fullfile(baseloc_pt, '*_Aff.csv'));

n_hc = length(FilesHC);
n_pt = length(FilesPT);

fprintf('Loading %d HC and %d Patient files...\n', n_hc, n_pt);

BigAffHC = zeros(214, 214, n_hc);
for s = 1:n_hc
    Aff = readtable(fullfile(baseloc_hc, FilesHC(s).name));
    BigAffHC(:,:,s) = table2array(Aff(:,2:end));
end

BigAffPT = zeros(214, 214, n_pt);
for s = 1:n_pt
    Aff = readtable(fullfile(baseloc_pt, FilesPT(s).name));
    BigAffPT(:,:,s) = table2array(Aff(:,2:end));
end

%% 3. Compute the network means
n_pairs = size(templateMat, 3);
BtwnNetworkConn = cell(n_pairs, 1);

for bb = 1:n_pairs

    mask = squeeze(templateMat(:,:,bb));
    tmpData = NaN(max(n_hc, n_pt), 2);

    for s = 1:n_hc
        subjMat = BigAffHC(:,:,s);
        tmpData(s,1) = mean(subjMat(mask == 1), 'omitnan');
    end

    for s = 1:n_pt
        subjMat = BigAffPT(:,:,s);
        tmpData(s,2) = mean(subjMat(mask == 1), 'omitnan');
    end

    BtwnNetworkConn{bb} = tmpData;
end

%% 4. Linear mixed effect model (setting the columns and calculation)
% Columns:
% 1  = p_raw
% 2  = FStat
% 3  = DF1
% 4  = DF2
% 5  = partial eta squared
% 6  = HC mean
% 7  = OUD mean
% 8  = HC SD
% 9  = OUD SD
% 10 = FDR p-value
% 11 = Model used: 1 = fitlme, 2 = fitlm fallback
% 12 = Hedges' g

BigRes = NaN(n_pairs, 12);

HC_ID = strcat("HC_", string(1:n_hc)');
PT_ID = strcat("OUD_", string(1:n_pt)');

for h = 1:n_pairs

    g1 = BtwnNetworkConn{h}(1:n_hc, 1);  % HC
    g2 = BtwnNetworkConn{h}(1:n_pt, 2);  % OUD

    BetweenNetworkEC = [g1; g2];

    % 0 = HC, 1 = OUD
    Group = [zeros(n_hc, 1); ones(n_pt, 1)];
    Subject_ID = [HC_ID; PT_ID];

    ECTable = table(BetweenNetworkEC, Group, Subject_ID);

    ECTable.Group = categorical(ECTable.Group);
    ECTable.Subject_ID = categorical(ECTable.Subject_ID);

    ECTable = rmmissing(ECTable);

    try
        lme = fitlme(ECTable, 'BetweenNetworkEC ~ Group + (1|Subject_ID)');
        lmeans = anova(lme);

        groupRow = strcmp(lmeans.Term, 'Group');

        p_raw = lmeans.pValue(groupRow);
        FStat = lmeans.FStat(groupRow);
        DF1 = lmeans.DF1(groupRow);
        DF2 = lmeans.DF2(groupRow);

        modelUsed = 1;

    catch ME
        warning('fitlme failed for pair %d: %s', h, ME.message);
        warning('Using fitlm fallback for this pair.');

        mdl = fitlm(ECTable, 'BetweenNetworkEC ~ Group');
        lmeans = anova(mdl);

        groupRow = strcmp(lmeans.Properties.RowNames, 'Group');

        p_raw = lmeans.pValue(groupRow);
        FStat = lmeans.F(groupRow);
        DF1 = lmeans.DF(groupRow);
        DF2 = mdl.DFE;

        modelUsed = 2;
    end

    %% Partial eta squared
    partial_eta2 = (FStat .* DF1) ./ ((FStat .* DF1) + DF2);

    %% Group descriptives
    HC_mean = mean(g1, 'omitnan');
    OUD_mean = mean(g2, 'omitnan');

    HC_sd = std(g1, 'omitnan');
    OUD_sd = std(g2, 'omitnan');

    %% Hedges' g
    n1 = sum(~isnan(g1));
    n2 = sum(~isnan(g2));

    pooled_sd = sqrt(((n1 - 1) * HC_sd^2 + (n2 - 1) * OUD_sd^2) / ...
                     (n1 + n2 - 2));

    if pooled_sd == 0 || isnan(pooled_sd)
        cohens_d = NaN;
        hedges_g = NaN;
    else
        % Positive g = OUD > HC
        % Negative g = HC > OUD
        cohens_d = (OUD_mean - HC_mean) / pooled_sd;

        % Small sample correction
        J = 1 - (3 / (4 * (n1 + n2) - 9));

        hedges_g = J * cohens_d;
    end

    %% The results will be stored in a BigRes table
    BigRes(h,1) = p_raw;
    BigRes(h,2) = FStat;
    BigRes(h,3) = DF1;
    BigRes(h,4) = DF2;
    BigRes(h,5) = partial_eta2;

    BigRes(h,6) = HC_mean;
    BigRes(h,7) = OUD_mean;
    BigRes(h,8) = HC_sd;
    BigRes(h,9) = OUD_sd;

    BigRes(h,11) = modelUsed;
    BigRes(h,12) = hedges_g;
end

%% 5. BH-FDR CORRECTION for multiple comparisons with significance at p<0.05.
m = n_pairs;
raw_p = BigRes(:,1);

[sorted_p, sort_idx] = sort(raw_p);

adj_p = sorted_p .* (m ./ (1:m)');

for i = m-1:-1:1
    adj_p(i) = min(adj_p(i), adj_p(i+1));
end

adj_p = min(adj_p, 1);

BigRes(sort_idx, 10) = adj_p;

%% 6. Print out to notify the results
fprintf('\n--- LMEM / LM BH-FDR Diagnostic for %s ---\n', taskType);
fprintf('%-5s %-30s %-12s %-12s %-10s %-10s %-10s %-10s\n', ...
    'Rank', 'Network Pair', 'Raw p', 'FDR p', 'FStat', 'pEta2', 'HedgesG', 'Model');
fprintf('%s\n', repmat('-', 1, 115));

for i = 1:m
    orig_idx = sort_idx(i);

    if BigRes(orig_idx,11) == 1
        modelLabel = 'fitlme';
    else
        modelLabel = 'fitlm';
    end

    fprintf('%-5d %-30s %-12.4f %-12.4f %-10.3f %-10.3f %-10.3f %-10s\n', ...
        i, ...
        TmpNetName{orig_idx}, ...
        sorted_p(i), ...
        adj_p(i), ...
        BigRes(orig_idx,2), ...
        BigRes(orig_idx,5), ...
        BigRes(orig_idx,12), ...
        modelLabel);
end

fprintf('\n');

%% 7. FINAL TABLE & SAVE
BigResTbl = table();

BigResTbl.NetworkPair = TmpNetName';
BigResTbl.p_raw = BigRes(:,1);
BigResTbl.p_fdr = BigRes(:,10);
BigResTbl.FStat = BigRes(:,2);
BigResTbl.DF1 = BigRes(:,3);
BigResTbl.DF2 = BigRes(:,4);
BigResTbl.Partial_Eta2 = BigRes(:,5);
BigResTbl.HedgesG = BigRes(:,12);
BigResTbl.HC_Mean = BigRes(:,6);
BigResTbl.HC_SD = BigRes(:,8);
BigResTbl.OUD_Mean = BigRes(:,7);
BigResTbl.OUD_SD = BigRes(:,9);

ModelUsed = strings(n_pairs,1);
ModelUsed(BigRes(:,11) == 1) = "fitlme";
ModelUsed(BigRes(:,11) == 2) = "fitlm_fallback";
BigResTbl.ModelUsed = ModelUsed;

sig_labels = repmat({''}, n_pairs, 1);

for i = 1:n_pairs
    if BigRes(i,10) < 0.05
        sig_labels{i} = '*';
    elseif BigRes(i,1) < 0.05 && BigRes(i,10) >= 0.05
        sig_labels{i} = 'trend';
    end
end

BigResTbl.Significance = sig_labels;

saveName = sprintf('Network_LMEM_Stats_%s.csv', taskType);
writetable(BigResTbl, fullfile(outputpath, saveName));

fprintf('Analysis for %s complete.\n', taskType);
fprintf('Results saved to: %s\n', fullfile(outputpath, saveName));
fprintf('Significant after FDR: %d / %d pairs\n', sum(BigRes(:,10) < 0.05), n_pairs);
fprintf('Trend-level raw p < .05, FDR p >= .05: %d / %d pairs\n', ...
    sum(BigRes(:,1) < 0.05 & BigRes(:,10) >= 0.05), n_pairs);
fprintf('fitlme models used: %d / %d\n', sum(BigRes(:,11) == 1), n_pairs);
fprintf('fitlm fallback models used: %d / %d\n', sum(BigRes(:,11) == 2), n_pairs);
