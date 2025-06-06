%Make notes about what the script is doing, number and comment like the lsNGC script

% Load BigResTbl and filter significant connections
BigResTbl = array2table(BigRes);
BigResTbl.Properties.VariableNames = {'r1', 'r2', 'p', 'FDRp', 'tstat', 'EffectSize'};

% Define the significance threshold
significance_threshold = 0.05;

% Filter the significant connections
significant_connections = BigResTbl(BigResTbl.FDRp < significance_threshold, :);

% Split into HC>HU and vice versa based on tstat
HCgBigResTbl = significant_connections(significant_connections.tstat > 0, :);
%HUgBigResTbl = significant_connections(significant_connections.tstat < 0, :);

% Specify regions to analyze
specified_regions = [201, 173, 174, 14, 110, 55];


% Display connections for specified regions in HC > HU
disp('Connections for specified regions in HC > HU:');
displayConnections(HCgBigResTbl, specified_regions);

% Display connections for specified regions in HU > HC
%disp('Connections for specified regions in HU > HC:');
%displayConnections(HUgBigResTbl, specified_regions);

% Function to compute connections
function connections = computeConnections(data, region)
    connections = [];
    for i = 1:height(data)
        if data.r1(i) == region || data.r2(i) == region
            if data.r1(i) == region
                target = data.r2(i);
            else
                target = data.r1(i);
            end
            connections = [connections; {region, target, data.EffectSize(i)}];
        end
    end
end

% Function to display connections for specified regions
function displayConnections(data, regions)
    for region = regions
        connections = computeConnections(data, region);
        if isempty(connections)
            disp(['No connections for region ', num2str(region)]);
        else
            total_connections = sum([connections{:, 3}]);
            for i = 1:size(connections, 1)
                connections{i, 3} = connections{i, 3} / total_connections; % Normalize
            end
            disp(['Connections for region ', num2str(region), ':']);
            disp(cell2table(connections, 'VariableNames', {'SourceRegion', 'TargetRegion', 'NormalizedCount'}));
        end
    end
end


