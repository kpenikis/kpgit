function [trainedClassifier, resubstitutionAccuracy] = CARTEST_trainClassifier
% trainClassifier(trainingData)
%  returns a trained classifier and its accuracy.
%  This code recreates the classification model trained in
%  Classification Learner app.
%
%   Input:
%       trainingData: the training data of same data type as imported
%        in the app (table or matrix).
%
%   Output:
%       trainedClassifier: a struct containing the trained classifier.
%        The struct contains various fields with information about the
%        trained classifier.
%
%       trainedClassifier.predictFcn: a function to make predictions
%        on new data. It takes an input of the same form as this training
%        code (table or matrix) and returns predictions for the response.
%        If you supply a matrix, include only the predictors columns (or
%        rows).
%
%       resubstitutionAccuracy: a double containing the accuracy in
%        percent. In the app, the History list displays this
%        overall accuracy score for each model.
%
%  Use the code to train the model with new data.
%  To retrain your classifier, call the function from the command line
%  with your original data or new data as the input argument trainingData.
%
%  For example, to retrain a classifier trained with the original data set
%  T, enter:
%    [trainedClassifier, resubstitutionAccuracy] = trainClassifier(T)
%
%  To make predictions with the returned 'trainedClassifier' on new data T,
%  use
%    yfit = trainedClassifier.predictFcn(T)
%
%  To automate training the same classifier with new data, or to learn how
%  to programmatically train classifiers, examine the generated code.

% Auto-generated by MATLAB on 02-Dec-2019 14:08:46


% Extract predictors and response
% This code processes the data into the right shape for training the
% classifier.


load carsmall
cartable = table(Acceleration, Cylinders, Displacement,Horsepower, Model_Year, MPG, Weight, Origin);

trainingData = cartable;

inputTable = trainingData;
predictorNames = {'Acceleration', 'Cylinders', 'Displacement', 'Horsepower', 'Model_Year', 'MPG', 'Weight'};
predictors = inputTable(:, predictorNames);
response = inputTable.Origin;
isCategoricalPredictor = [false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true);
classificationSVM = fitcecoc(...
    predictors, ...
    response, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'ClassNames', ['F' 'r' 'a' 'n' 'c' 'e' ' '; 'G' 'e' 'r' 'm' 'a' 'n' 'y'; 'I' 't' 'a' 'l' 'y' ' ' ' '; 'J' 'a' 'p' 'a' 'n' ' ' ' '; 'S' 'w' 'e' 'd' 'e' 'n' ' '; 'U' 'S' 'A' ' ' ' ' ' ' ' ']);

% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
svmPredictFcn = @(x) predict(classificationSVM, x);
trainedClassifier.predictFcn = @(x) svmPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'Acceleration', 'Cylinders', 'Displacement', 'Horsepower', 'Model_Year', 'MPG', 'Weight'};
trainedClassifier.ClassificationSVM = classificationSVM;
trainedClassifier.About = 'This struct is a trained classifier exported from Classification Learner R2016b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedClassifier''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Extract predictors and response
% This code processes the data into the right shape for training the
% classifier.
inputTable = trainingData;
predictorNames = {'Acceleration', 'Cylinders', 'Displacement', 'Horsepower', 'Model_Year', 'MPG', 'Weight'};
predictors = inputTable(:, predictorNames);
response = inputTable.Origin;
isCategoricalPredictor = [false, false, false, false, false, false, false];

% Compute resubstitution accuracy
resubstitutionAccuracy = 1 - resubLoss(trainedClassifier.ClassificationSVM, 'LossFun', 'ClassifError');

% Compute resubstitution predictions and scores
[resubstitutionPredictions, resubstitutionScores] = predict(trainedClassifier.ClassificationSVM, predictors);
