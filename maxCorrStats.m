function [max_r, min_p, r_95th, adj_p, rval_null, pval_null] = ...
    maxCorrStats(x, y, z, test_r, test_p, type, rows, iterations)
% This function performs partial correlations between x and a specified
% number of random permutations of y, controlling for the variables in z.
% It creates an empirical null distribution against which actual test
% statistics can be compared. It returns the variables in test_r and test_p
% that are deemed significant according to the empirical null distribution.
% INPUT:
% x          = array of size nparticipants * nvariables containing data 
%              from outcome measures to be correlated with data in y.
% y          = array of size nparticipants * 1 containing data to be
%              correlated with outcome measures in x.
% z          = array of size nparticipants * nvariables containing data 
%              from variables used to control/adjust the correlations 
%              between x and y.
% test_r     = table containing r/rho values obtained from partial
%              correlations of x and y, controlling for z.
% test_p     = table containing p values obtained from partial correlations
%              of x and y, controlling for z.
% type       = string indicating type of correlation to carry out using 
%              partialcorr.m. Options: 'Pearson' to carry out Pearson
%              partial correlation; 'Spearman' to carry out Spearman 
%              partial correlation. *NOTE: The type of partial correlation
%              selected here must be the same partial correlation as
%              carried out to obtain the values in test_r and test_p
% rows       = string indicating how to deal with missing values. Options:
%              'all' to use all rows regardless of NaNs; 'complete' to use
%              listwise deletion; 'pairwise' to use pairwise deletion
% iterations = double indicating how many times to perform partial
%              correlations with randomly shuffled y
% OUTPUT:
% max_rho    = table containing variable names and r/rho values of 
%              variables with r/rho values that are greater than r/rho 
%              value in 95th percentile in randomly generated null
%              distribution of r/rho values.
% min_p      = table containing variable names and p values of variables
%              with p values that are smaller than p value in 5th 
%              percentile in randomly generated null distribution of p
%              values.
% r_95th     = 95th percentile of r/rho values (threshold above which r/rho
%              values are considered statistically significant
% adj_p      = adjusted pvalue accounting for multiple comparisons and
%              intercorrelations of variables in x. Pvalues below adj_p can
%              be considered statistically significant.
% rval_null  = nvariables * niterations array containing r/rho values of
%              correlations between x and random permutations of y
% pval_null  = nvariables * niterations array containing p values of
%              correlations between x and random permutations of y

% Author: Rory Boyle rorytboyle@gmail.com
% Date: 29/09/18

%% 1) Parse and clean input
% check x, y, and z are numerical arrays
if ~isnumeric(x) | ~isnumeric(y) | ~isnumeric(z)
    error('x and/or y and/or z are not provided as numerical arrays')
end

% check test_r and test_p are tables
if ~istable(test_r) | ~istable(test_p)
    error('test_r and/or test_p are not provided as matlab tables, please structure these vars as tables')
end

% check type is provided correctly
if ~strcmp(type, 'Pearson') & ~strcmp(type, 'Spearman')
    error("type is not 'Pearson' or 'Spearman', these are the only two options - make sure type is spelt correctly")
end

% check rows are provided correctly
if ~strcmp(rows, 'all') & ~strcmp(rows, 'complete') & ...
        ~strcmp(rows, 'pairwise')
    error("rows is not 'all' or 'complete' or 'pairwise', these are the only three options - make sure rows is spelt correctly")
end

% check iterations is a number
if ~strcmp(class(iterations), 'double')
    error('iterations is not provided as a double')
end
%% 2) Generate null distributions 
% init pval and rho matrix 
pval_null = nan(size(x,2), iterations);

rval_null = nan(size(x,2), iterations);% rval refers to r if type = Pearson
                                       % and rho if type = Spearman
                                       
% randomly shuffle y and perform partial correlations (repeat for the 
% number of times specified by the input 'iterations')                                        
for i = 1:iterations
    % shuffle y
    y_shuffled = y(randperm(length(y)));
    
    % do partial correlation
    [rval_null(:,i), pval_null(:,i)] = partialcorr(x, y_shuffled,...
        z, 'Type', type, 'Rows', rows);
end

%% 3) Calculate percentile values

% first use absolute value of rho to account for both positive and negative
% correlations
rval_null = abs(rval_null);

% get vector of max r/rho values and min p values %%% Added by Rob Whelan
maxRhoVals = max(rval_null);
minPvals = min(pval_null);

% get percentiles
r_95th = prctile(maxRhoVals, 95); % max r/rho
adj_p = prctile(minPvals, 5);   % min p

%% 4)compare to actual test statistics
% return variables and their rvalues for those with rvalues greater than/ 
% equal to value in 95th percentile of null model of r values 
% use absolute values of test r/rho to account for both positive and
% negative correlations
max_r = test_r(abs(table2array(test_r))>=r_95th,:);

% return variables and their pvalues for those with pvalues smaller than/ 
% equal to value in 5th percentile of null model of p values
min_p = test_p(table2array(test_p)<=adj_p,:);

end
