function [Parameters] = learnCIMeasureParams()
% Parameters for learnCIMeasure_noisyor()
% OUTPUT 
%  Parameters - struct - The struct contains the following fields:
%                % These parameters are user defined (can be modified)
%                    1. nPop: Size of population
%                    2. sigma: Sigma of Gaussians in fitness function
%                    3. nIterations: Number of iterations
%                    4. eta:Percentage of time to make small-scale mutation
%                    5. sampleVar: Variance around sample mean
%                    6. mean: mean of ci in fitness function. Is always set to 1 if the positive label is "1".
%                    7. analysis: if ="1", record all intermediate results
%
               
%Parameters.nPop = 30; %Size of population (parallel processing not implemented)
Parameters.eta = 0.5; %Percentage of time to make small-scale mutation
Parameters.analysis = 1; % if ="1", record all intermediate results
Parameters.exaustiveSearchThresh = 500; % count threshold for number of repeated samples
Parameters.fitnessUpdateThresh = 100; % count threshold for number of times new BFM samples do not improve over past iterations
end

 