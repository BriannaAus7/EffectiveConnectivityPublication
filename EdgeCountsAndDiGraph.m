% This script is for the dominant network of influence analysis. 
% The dominant network of influence as a measure of hierarchal control was based on the highest connectivity counts from the edge patterns. 
% The highest connectivity count indicates a networks large influence across other networks. 
% To determine the dominant network, the frequency of each network’s appearance was counted across all edge patterns, excluding self-connections (e.g., controlcontrol). 
% Connections were split into HC>MD and MD>HC based on z-statistic direction.

% 1. Load BigResTbl and filter significant connections
BigResTbl = array2table(BigRes);
BigResTbl.Properties.VariableNames = {'r1', 'r2', 'p', 'FDRp', 'zstat', 'EffectSize'};

% 2. Define the significance threshold
significance_threshold = 0.05;

% 3. Filter the significant connections
significant_connections = BigResTbl(BigResTbl.FDRp < significance_threshold, :);

% 4. Split the data into HC > MD and MD > HC based on zstat
HC_greater_MD = significant_connections(significant_connections.zstat > 0, :);
MD_greater_HC = significant_connections(significant_connections.zstat < 0, :);

% 5. Call the function with the HC > MD data
[WorkingTable_HC, PropRedNetList_HC, RedNetNames_HC, RegionCounts_HC, EdgePatterns_HC] = ComputeNetworkDegree(HC_greater_MD);

% 6. Call the function with the MD > HC data
[WorkingTable_MD, PropRedNetList_MD, RedNetNames_MD, RegionCounts_MD, EdgePatterns_MD] = ComputeNetworkDegree(MD_greater_HC);

% 7. Convert RegionCounts and EdgePatterns maps to tables for better readability (HC > MD)
RegionCountsTable_HC = cell2table([keys(RegionCounts_HC)' values(RegionCounts_HC)'], 'VariableNames', {'Region', 'Count'});
EdgePatternsTable_HC = cell2table([keys(EdgePatterns_HC)' values(EdgePatterns_HC)'], 'VariableNames', {'EdgePattern', 'Count'});

% 8. Sort the tables by count in descending order (HC > MD)
RegionCountsTable_HC = sortrows(RegionCountsTable_HC, 'Count', 'descend');
EdgePatternsTable_HC = sortrows(EdgePatternsTable_HC, 'Count', 'descend');

% 9. Display the results (HC > MD)
disp('Region Counts (HC > MD, sorted by count):');
disp(RegionCountsTable_HC);
disp('Edge Patterns (HC > MD, sorted by count):');
disp(EdgePatternsTable_HC);

% 10. Convert RegionCounts and EdgePatterns maps to tables for better readability (MD > HC)
RegionCountsTable_MD = cell2table([keys(RegionCounts_MD)' values(RegionCounts_MD)'], 'VariableNames', {'Region', 'Count'});
EdgePatternsTable_MD = cell2table([keys(EdgePatterns_MD)' values(EdgePatterns_MD)'], 'VariableNames', {'EdgePattern', 'Count'});

% 11. Sort the tables by count in descending order (MD > HC)
RegionCountsTable_MD = sortrows(RegionCountsTable_MD, 'Count', 'descend');
EdgePatternsTable_MD = sortrows(EdgePatternsTable_MD, 'Count', 'descend');

% 12. Display the results (MD > HC)
disp('Region Counts (MD > HC, sorted by count):');
disp(RegionCountsTable_MD);
disp('Edge Patterns (MD > HC, sorted by count):');
disp(EdgePatternsTable_MD);

% 13. Define the colormap for each network
network_labels = {'Vis', 'SoMat', 'DorsAttn', 'SalVent', 'Limb', 'Cntrl', 'DMN', 'TempPar', 'SubCor', 'VMN'};
num_networks = length(network_labels);

% 14. Define colors for each network
colors = [
    0.5 0 0.5;    % Vis - Purple
    1 0.75 0.8;   % SoMat - Pink
    0 1 0;        % DorsAttn - Green
    0.6 0.3 0.2;  % SalVent - Brown
    1 0.65 0;     % Limb - Orange
    1 0 0;        % Cntrl - Red
    0 0 0.5;      % DMN - Dark Blue
    0 0.5 0.5;    % TempPar - Teal
    0 0 0;        % SubCor - Black
    1 1 0         % VMN - Yellow
];

network_colors = containers.Map(network_labels, num2cell(colors, 2));

% 15. Normalize weights function
normalize_weights = @(weights) (weights - min(weights)) / (max(weights) - min(weights));

% 16. Prepare data for graph visualization (HC > MD)
connection_matrix_HC = zeros(num_networks);

% 17. Corrected loop for HC > MD
keys_EdgePatterns_HC = keys(EdgePatterns_HC);
values_EdgePatterns_HC = values(EdgePatterns_HC);

for i = 1:length(keys_EdgePatterns_HC)
    pattern = strsplit(keys_EdgePatterns_HC{i}, '-');
    pattern = str2double(pattern);
    count = values_EdgePatterns_HC{i};
    connection_matrix_HC(pattern(1), pattern(2)) = count;
end

% 18. Create directed graph (HC > MD)
G_HC = digraph(connection_matrix_HC, network_labels);

% 19. Plot the graph (HC > MD) with consistent colors and no edge labels
plot_graph(G_HC, network_labels, network_colors, 'Network Effective Connectivity (HC > MD)');

% 20. Prepare data for graph visualization (MD > HC)
connection_matrix_MD = zeros(num_networks);

% 21. Corrected loop for MD > HC
keys_EdgePatterns_MD = keys(EdgePatterns_MD);
values_EdgePatterns_MD = values(EdgePatterns_MD);

for i = 1:length(keys_EdgePatterns_MD)
    pattern = strsplit(keys_EdgePatterns_MD{i}, '-');
    pattern = str2double(pattern);
    count = values_EdgePatterns_MD{i};
    connection_matrix_MD(pattern(1), pattern(2)) = count;
end

% 22. Create directed graph (MD > HC)
G_MD = digraph(connection_matrix_MD, network_labels);

% 23. Plot the graph (MD > HC) with consistent colors and no edge labels
plot_graph(G_MD, network_labels, network_colors, 'Network Effective Connectivity (MD > HC)');



function plot_graph(G, network_labels, network_colors, title_str)
    % Plot the graph
    figure;
    % Remove 'EdgeLabel' to hide connection numbers
    p = plot(G, 'Layout', 'circle');
    title(title_str);
    xlabel('Networks');
    ylabel('Connection Strength');
    
    % 24. Customize edge thickness based on weights, using non-linear scaling for contrast
    weights = G.Edges.Weight;
    normalized_weights = (weights - min(weights)) / (max(weights) - min(weights)); % Normalize weights to [0, 1]
    adjusted_weights = normalized_weights .^ 2; % Apply non-linear scaling to enhance contrast
    p.LineWidth = 1 + 8 * adjusted_weights; % Increase range for more visual distinction
    p.EdgeColor = [1 0.5 0]; % Set all edges to red

    % 25. Adjust font size for node labels
    p.NodeFontSize = 14; % Adjust this value to make node labels larger
    
    p.ArrowSize = 12; % Adjust this value to make arrows larger (default is usually 7)

    % 26. Apply colors to nodes only
    for i = 1:length(network_labels)
        highlight(p, i, 'NodeColor', network_colors(network_labels{i}));
    end
end

function [WorkingTable, PropRedNetList, RedNetNames, RegionCounts, EdgePatterns] = ComputeNetworkDegree(BigResTbl)
    % Load network template data
    load TemplateNets.mat

    % 27. Initialize table to store results
    WorkingTable = table();
    WorkingTable.r1 = BigResTbl.r1;
    WorkingTable.r2 = BigResTbl.r2;
    WorkingTable.Network1 = zeros(size(BigResTbl, 1), 1);
    WorkingTable.Network2 = zeros(size(BigResTbl, 1), 1);

    % 28. Initialize counts and patterns
    VisCount = 0; SoMatCount = 0; DorsAttnCount = 0; SalVentCount = 0;
    LimbCount = 0; ContCount = 0; DMNCount = 0; TmpParCount = 0; SubCorCount = 0; VMNCount = 0;
    RegionCounts = containers.Map('KeyType', 'double', 'ValueType', 'double');
    EdgePatterns = containers.Map('KeyType', 'char', 'ValueType', 'double');

    % 29. Counting ROI1
    for ii = 1:size(WorkingTable, 1)
        roi1 = WorkingTable.r1(ii);

        if isKey(RegionCounts, roi1)
            RegionCounts(roi1) = RegionCounts(roi1) + 1;
        else
            RegionCounts(roi1) = 1;
        end

        if ismember(roi1, Vis)
            VisCount = VisCount + 1;
            WorkingTable.Network1(ii) = 1;
        elseif ismember(roi1, SoMat)
            SoMatCount = SoMatCount + 1;
            WorkingTable.Network1(ii) = 2;
        elseif ismember(roi1, DorsAttn)
            DorsAttnCount = DorsAttnCount + 1;
            WorkingTable.Network1(ii) = 3;
        elseif ismember(roi1, SalVent)
            SalVentCount = SalVentCount + 1;
            WorkingTable.Network1(ii) = 4;
        elseif ismember(roi1, Limb)
            LimbCount = LimbCount + 1;
            WorkingTable.Network1(ii) = 5;
        elseif ismember(roi1, Cont)
            ContCount = ContCount + 1;
            WorkingTable.Network1(ii) = 6;
        elseif ismember(roi1, DMN)
            DMNCount = DMNCount + 1;
            WorkingTable.Network1(ii) = 7;
        elseif ismember(roi1, TmpPar)
            TmpParCount = TmpParCount + 1;
            WorkingTable.Network1(ii) = 8;
        elseif ismember(roi1, SubCor)
            SubCorCount = SubCorCount + 1;
            WorkingTable.Network1(ii) = 9;
        elseif ismember(roi1, VMN)
            VMNCount = VMNCount + 1;
            WorkingTable.Network1(ii) = 10;
        else
            disp(strcat('No match for ROI1: ', num2str(roi1)));
        end
    end

    % 30. Counting ROI2
    for ii = 1:size(WorkingTable, 1)
        roi2 = WorkingTable.r2(ii);

        if isKey(RegionCounts, roi2)
            RegionCounts(roi2) = RegionCounts(roi2) + 1;
        else
            RegionCounts(roi2) = 1;
        end

        if ismember(roi2, Vis)
            VisCount = VisCount + 1;
            WorkingTable.Network2(ii) = 1;
        elseif ismember(roi2, SoMat)
            SoMatCount = SoMatCount + 1;
            WorkingTable.Network2(ii) = 2;
        elseif ismember(roi2, DorsAttn)
            DorsAttnCount = DorsAttnCount + 1;
            WorkingTable.Network2(ii) = 3;
        elseif ismember(roi2, SalVent)
            SalVentCount = SalVentCount + 1;
            WorkingTable.Network2(ii) = 4;
        elseif ismember(roi2, Limb)
            LimbCount = LimbCount + 1;
            WorkingTable.Network2(ii) = 5;
        elseif ismember(roi2, Cont)
            ContCount = ContCount + 1;
            WorkingTable.Network2(ii) = 6;
        elseif ismember(roi2, DMN)
            DMNCount = DMNCount + 1;
            WorkingTable.Network2(ii) = 7;
        elseif ismember(roi2, TmpPar)
            TmpParCount = TmpParCount + 1;
            WorkingTable.Network2(ii) = 8;
        elseif ismember(roi2, SubCor)
            SubCorCount = SubCorCount + 1;
            WorkingTable.Network2(ii) = 9;
        elseif ismember(roi2, VMN)
            VMNCount = VMNCount + 1;
            WorkingTable.Network2(ii) = 10;
        else
            disp(strcat('No match for ROI2: ', num2str(roi2)));
        end

        % 31. Track edge patterns
        edge_pattern = strcat(num2str(WorkingTable.Network1(ii)), '-', num2str(WorkingTable.Network2(ii)));
        if isKey(EdgePatterns, edge_pattern)
            EdgePatterns(edge_pattern) = EdgePatterns(edge_pattern) + 1;
        else
            EdgePatterns(edge_pattern) = 1;
        end
    end

    % 32. Compute degree per 10 networks.
    RedNetList = [VisCount, SoMatCount, DorsAttnCount, SalVentCount, LimbCount, ContCount, DMNCount, TmpParCount, SubCorCount, VMNCount]';
    RedNetNames = {"Vis", "SoMat", "DorsAttn", "SalVent", "Limb", "Cntrl", "DMN", "TempPar", "SubCor", "VMN"}';

    PropRedNetList = RedNetList / sum(RedNetList);
end
