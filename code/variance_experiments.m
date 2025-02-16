%% Main simulation script for Gaussian or Gamma (shifted to median or mean = 0)
clear; close all;

%% ===== User Settings =====
% Choose distribution type: 'gaussian' or 'gamma'
distType = 'gamma';  % Change to 'gaussian' for Gaussian simulation
% Choose centering method for Gamma: 'median' or 'mean'
centerType = 'mean';  % Set to 'mean' if you want to shift by the mean
nExperiments = 500;  % Number of experiments per parameter value  
n = 50;               % Sample size per group
nPermutations = 100; % Number of permutations in the permutation test

% For Gaussian, we loop over sigma (standard deviation) values.
% For Gamma, we loop over scale (θ) values with a fixed shape.
if strcmpi(distType, 'gaussian')
    paramValues = 0.2:0.5:1.5;  % These are σ values.
else
    % For Gamma: fix a shape parameter and loop over scale values.
    fixedShape = 2;            % Fixed shape parameter (a)
    paramValues = 0.2:0.5:4.0;  % These will be used as the scale (θ) values.
end

% skewness of 2 top right plot
% compute averages in the 2d plot
% plot gamma densities and p value
% increase the number of permutations with the skewness
% plot two seperate histograms for group 2 data at the bottom

% Settings for p-value histogram: bins of width 0.01 from 0 to 1
edges = 0:0.01:1;
binCenters = edges(1:end-1) + diff(edges)/2;

% Colors for plotting (one per parameter value):
colors = lines(length(paramValues));

%% Create a figure with six subplots:
figure;

%% --- Subplot 1: Histogram of p-values ---
subplot(3,2,1); % Create a 3x2 subplot grid, use the first subplot
hold on;
legendEntries = strings(length(paramValues), 1);
allPValues = {}; % Store p-values for each parameter value
allGroup2Data = {}; % Store Group 2 data for each parameter value

for idx = 1:length(paramValues)
    paramVal = paramValues(idx);
    pvalues = zeros(nExperiments, 1);
    group2Data = cell(nExperiments, 1); % Store Group 2 data for each experiment
    
    for iExp = 1:nExperiments
        % Group 1: all zeros.
        group1 = zeros(n, 1);
        
        % Group 2: generate data according to the selected distribution.
        if strcmpi(distType, 'gaussian')
            % For Gaussian: generate with mean = 0 and standard deviation = paramVal.
            group2 = normrnd(0, paramVal, [n, 1]);
        else
            % For Gamma: generate gamma data with fixed shape and scale = paramVal.
            % Shift by subtracting the chosen center so that the data have center = 0.
            if strcmpi(centerType, 'median')
                centerVal = gaminv(0.5, fixedShape, paramVal);
            else  % use mean
                centerVal = fixedShape * paramVal;
            end
            group2 = gamrnd(fixedShape, paramVal, [n, 1]) - centerVal;
        end
        
        % Store Group 2 data
        group2Data{iExp} = group2;
        
        % Compute the two-tailed p-value using a permutation test (difference of means).
        p = permutationMeanTest(group1, group2, nPermutations);
        pvalues(iExp) = p;
    end
    
    % Store p-values and Group 2 data for later use
    allPValues{idx} = pvalues;
    allGroup2Data{idx} = group2Data;
    
    % Compute histogram counts for the p-values using the specified bins.
    counts = histcounts(pvalues, edges);
    
    % Plot the histogram as a bar plot.
    bar(binCenters, counts, 'FaceColor', colors(idx, :), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.6);
    
    % Create an appropriate legend entry.
    if strcmpi(distType, 'gaussian')
        legendEntries(idx) = sprintf('\\sigma = %.1f', paramVal);
    else
        if strcmpi(centerType, 'median')
            legendEntries(idx) = sprintf('Scale = %.1f, Shape = %.1f (Median)', paramVal, fixedShape);
        else
            legendEntries(idx) = sprintf('Scale = %.1f, Shape = %.1f (Mean)', paramVal, fixedShape);
        end
    end
end
xlabel('p-value');
ylabel('Frequency');
if strcmpi(distType, 'gaussian')
    title('Histogram of p-values (Gaussian)');
else
    if strcmpi(centerType, 'median')
        title('Histogram of p-values (Gamma, Median)');
    else
        title('Histogram of p-values (Gamma, Mean)');
    end
end
legend(legendEntries, 'Location', 'northwest');
grid on;
hold off;

%% --- Subplot 2: Theoretical Density Functions ---
subplot(3,2,2); % Use the second subplot
hold on;
legendEntries2 = strings(length(paramValues), 1);
for idx = 1:length(paramValues)
    paramVal = paramValues(idx);
    if strcmpi(distType, 'gaussian')
        % For Gaussian: mean = 0, standard deviation = paramVal.
        sigma = paramVal;
        x = linspace(-4*sigma, 4*sigma, 1000);
        y = normpdf(x, 0, sigma);
        legendEntries2(idx) = sprintf('\\sigma = %.1f', sigma);
    else
        % For Gamma: fixed shape and scale = paramVal.
        theta = paramVal;
        if strcmpi(centerType, 'median')
            centerVal = gaminv(0.5, fixedShape, theta);
        else
            centerVal = fixedShape * theta;
        end
        sigma_gamma = sqrt(fixedShape) * theta;  % standard deviation
        
        % Define an x-range: starting at -centerVal (shifted support) and going up
        % To capture most of the mass, use the 99th percentile of the unshifted gamma.
        xMin = -centerVal;
        xMax = gaminv(0.99, fixedShape, theta) - centerVal;
        x = linspace(xMin, xMax, 1000);
        % Evaluate the density of the shifted gamma:
        y = gampdf(x + centerVal, fixedShape, theta);
        if strcmpi(centerType, 'median')
            legendEntries2(idx) = sprintf('Scale = %.1f, Shape = %.1f (Median)', theta, fixedShape);
        else
            legendEntries2(idx) = sprintf('Scale = %.1f, Shape = %.1f (Mean)', theta, fixedShape);
        end
    end
    plot(x, y, 'Color', colors(idx, :), 'LineWidth', 2);
end
xlabel('x');
ylabel('Probability Density');
if strcmpi(distType, 'gaussian')
    title('Gaussian Densities');
else
    if strcmpi(centerType, 'median')
        title('Gamma Densities (Median)');
    else
        title('Gamma Densities (Mean)');
    end
end
legend(legendEntries2, 'Location', 'northeast');
grid on;
hold off;

%% --- Subplot 3: Surface Plot: Type I Error Inflation for Gamma (shifted accordingly) ---
subplot(3,2,3); % Use the third subplot
% Define a grid for shape and scale parameters:
shapeVec = 1:0.5:10;         % e.g., shape values from 1 to 10
scaleVec = 0.5:0.5:4.0;      % e.g., scale values from 0.5 to 4.0
nExpSurface = 500;  % Number of experiments per grid cell (reduce if needed)

% Preallocate matrices for inflation, standard deviation, and skewness.
inflationMat = zeros(length(shapeVec), length(scaleVec));
stdMat = zeros(length(shapeVec), length(scaleVec));
skewMat = zeros(length(shapeVec), length(scaleVec));

for i = 1:length(shapeVec)
    for j = 1:length(scaleVec)
        shapeVal = shapeVec(i);
        scaleVal = scaleVec(j);
        if strcmpi(centerType, 'median')
            centerVal = gaminv(0.5, shapeVal, scaleVal);
        else
            centerVal = shapeVal * scaleVal;
        end
        % Compute theoretical standard deviation and skewness.
        sigmaVal = sqrt(shapeVal) * scaleVal;
        skewVal = 2 / sqrt(shapeVal);
        
        stdMat(i, j) = sigmaVal;
        skewMat(i, j) = skewVal;
        
        % Simulate experiments for this parameter combination.
        pvals = zeros(nExpSurface, 1);
        for k = 1:nExpSurface
            group1 = zeros(n, 1);
            % Shift Gamma data by the chosen center.
            group2 = gamrnd(shapeVal, scaleVal, [n, 1]) - centerVal;
            pvals(k) = permutationMeanTest(group1, group2, nPermutations);
        end
        % Inflation: fraction of experiments with p < 0.05
        inflationMat(i, j) = mean(pvals < 0.05);
    end
end

surf(stdMat, skewMat, inflationMat);
xlabel('Standard Deviation');
ylabel('Skewness');
zlabel('Type I Error Inflation');
if strcmpi(centerType, 'median')
    title('Type I Error Inflation (Median)');
else
    title('Type I Error Inflation (Mean)');
end
colorbar;
shading interp;

%% --- Subplot 4: Skewness vs Type I Error Rate ---
subplot(3,2,4); % Use the fourth subplot
scatter(skewMat(:), inflationMat(:), 50, 'filled');
xlabel('Skewness');
ylabel('Type I Error Inflation');
title('Skewness vs. Type I Error Rate');
grid on;

%% --- Subplot 5: Data from Group 2 Corresponding to Specific p-value Ranges ---
subplot(3,2,[5,6]); % Use the fifth and sixth subplots for this plot
hold on;

% Define p-value ranges
pValueRanges = [0, 0.01, 0.05, 0.25, 0.5, 0.75, 1];
rangeLabels = {'0-1%', '1-5%', '5-25%', '25-50%', '50-75%', '75-100%'};

% Preallocate data containers for each range
group2DataByRange = cell(1, length(rangeLabels));

% Loop through all stored p-values and Group 2 data
for idx = 1:length(paramValues)
    pvalues = allPValues{idx};
    group2Data = allGroup2Data{idx};

    for iExp = 1:nExperiments
        p = pvalues(iExp);
        group2 = group2Data{iExp};

        % Assign Group 2 data to the appropriate range
        for r = 1:length(rangeLabels)
            if p >= pValueRanges(r) && p < pValueRanges(r+1)
                group2DataByRange{r} = [group2DataByRange{r}; group2];
                break;
            end
        end
    end
end

% Combine all data to determine the overall range
allData = [];
for r = 1:length(rangeLabels)
    allData = [allData; group2DataByRange{r}];
end

% Determine the bin edges or bin width
minValue = min(allData);
maxValue = max(allData);
numBins = 100; % Specify the number of bins you want
binEdges = linspace(minValue, maxValue, numBins + 1);

% Plot data for each range
for r = 1:length(rangeLabels)
    if ~isempty(group2DataByRange{r})
        histogram(group2DataByRange{r}, binEdges, 'Normalization', 'pdf', ...
                  'FaceColor', colors(r,:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0.6);
        hold on;
    end
end

xlabel('Group 2 Data Values');
ylabel('Density');
title('Group 2 Data by p-value Range');
legend(rangeLabels, 'Location', 'northeast');
xlim([-12, 12])
grid on;
hold off;

%% --- Local Function: Permutation Test Using Difference of Means ---
function p = permutationMeanTest(group1, group2, nPermutations)
    % Calculate the observed difference of means.
    obsDiff = mean(group1) - mean(group2);
    
    % Combine both groups.
    combined = [group1; group2];
    n = length(group1);  % assumes equal group sizes
    
    % Preallocate an array for permutation differences.
    permDiff = zeros(nPermutations, 1);
    
    for j = 1:nPermutations
        % Randomly permute the combined data.
        permutedIdx = randperm(length(combined));
        
        % Divide the permuted data into two groups.
        permGroup1 = combined(permutedIdx(1:n));
        permGroup2 = combined(permutedIdx(n+1:end));
        
        % Compute the difference of means for this permutation.
        permDiff(j) = mean(permGroup1) - mean(permGroup2);
    end
    
    % The two-tailed p-value is the fraction of permuted differences (in absolute value)
    % that are at least as large as the observed absolute difference.
    p = mean(abs(permDiff) >= abs(obsDiff));
end