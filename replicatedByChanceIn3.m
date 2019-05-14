function [probSigIn3, timesSigIn3] = replicatedByChanceIn3(x_1, x_2, x_3, y_1, y_2,...
    y_3, z_1, z_2, z_3, type, rows, iterations, p_1, p_2, p_3)
% This function performs partial correlations between x and a specified
% number of random permutations of y, controlling for the variables in z.
% It repeats these correlations a number of times specified by iterations 
% and returns a p-value representing the probability of the findings being 
% replicated across 3 datasets by chance (i.e. the number of times the 
% random/null model p-values were smaller than the actual pvalues in all 3 
% datasets, divided by the number of iterations)
%
% INPUT:
% x_1        = column vector of data from dataset 1 to be correlated with
%              data in y_1.
% x_2        = column vector of data from dataset 2 to be correlated with
%              data in y_2.
% x_3        = column vector of data from dataset 3 to be correlated with
%              data in y_3.
% y_1        = column vector of data from dataset 1 to be correlated with
%              data in x_1
% y_2        = column vector of data from dataset 2 to be correlated with
%              data in x_2
% y_3        = column vector of data from dataset 2 to be correlated with
%              data in x_3
% z_1        = array of size nparticipants * nvariables containing data 
%              from variables in dataset 1 used to control/adjust the
%              correlations between x_1 and y_1.
% z_2        = array of size nparticipants * nvariables containing data 
%              from variables in dataset 2 used to control/adjust the
%              correlations between x_2 and y_2.
% z_3        = array of size nparticipants * nvariables containing data 
%              from variables in dataset 1 used to control/adjust the
%              correlations between x_3 and y_3.
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
% p_1        = pvalue of actual partial correlation between x_1 and y_1
%              controlling for z_1
% p_2        = pvalue of actual partial correlation between x_2 and y_2
%              controlling for z_2
% p_3        = pvalue of actual partial correlation between x_3 and y_3
%              controlling for z_3
% OUTPUT:
% timesSigIn3 =double containing the number of iterations in which the
%              pvalues from random partial correlations were < actual
%              pvalues from real data in all 3 datasets
% probSigIn3 = double representing the probability of a finding being
%              replicated across 3 datasets by change. it is equivalent to 
%              timesSigIn3 divided by iterations

% Authors: Rory Boyle & Rob Whelan - rorytboyle@gmail.com
% Date: 10/10/18
% Update: 01/11/18 - Edited code to only look across all 3 datasets.
% Previous version replicatedByChance_old looked at likelihoods in 2/3 and
% 3/3 datasets. Also edited to calculate number of iterations where all
% null p values were less than all actual p values (regardless of dataset).
% Update: 28/03/2019 - Changed code to calculate probability of replicating
% significant finding in 3 datasets instead of percentage.

%% 1) Parse and clean input
% check x, y, and z are numerical arrays
if ~isnumeric(x_1) | ~isnumeric(y_1) | ~isnumeric(z_1) | ...
        ~isnumeric(x_2) | ~isnumeric(y_2) | ~isnumeric(z_2) | ...
        ~isnumeric(x_3) | ~isnumeric(y_3) | ~isnumeric(z_3)
    error('x and/or y and/or z are not provided as numerical arrays')
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

%% 2) Perform partial correlations of shuffled y values with x (controlling for z)
% init var to store p values - column 1 = dataset 1, column 2 = dataset 2, 
% column 3 = dataset 3
null_pvals = zeros(iterations,3);

% randomly shuffle y_1, y_2, and y_3 and perform partial correlations (repeat for
% the number of times specified by the input 'iterations')                                        
for i = 1:iterations

    % shuffle y values
    y_1_shuffled = y_1(randperm(length(y_1)));
    y_2_shuffled = y_2(randperm(length(y_2)));
    y_3_shuffled = y_3(randperm(length(y_3)));
    
    % do partial correlation dataset 1
    [~, null_pvals(i,1)] = partialcorr(x_1, y_1_shuffled,...
        z_1, 'Type', type, 'Rows', rows);
    
    % do partial correlation dataset 2
    [~, null_pvals(i,2)] = partialcorr(x_2, y_2_shuffled,...
        z_2, 'Type', type, 'Rows', rows);
   
    % do partial correlation dataset 3
    [~, null_pvals(i,3)] = partialcorr(x_3, y_3_shuffled,...
        z_3, 'Type', type, 'Rows', rows);
end

%% 3) Calculate probability of replicating by chance

% sort actual p values
actual_pvals = sort([p_1, p_2, p_3]);

% sort null pvalues for each iteration (i.e. sort each row)
sort_null_pvals=sort(null_pvals,2);

% count number of iterations in which all null (random) p values are
% smaller than the actual p values
timesSigIn3=sum(all(sort_null_pvals<actual_pvals,2));

% calculate and return p-value representing probability of replicating 
% findings by chance
probSigIn3 = timesSigIn3/iterations;
end
