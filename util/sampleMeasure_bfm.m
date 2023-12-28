function [measure] = sampleMeasure_bfm(eta, nSources, lowerindex, upperindex,  prevSample, IndexsubsetToChange)
%sampleMeasure - sampling a new measure
% If the number of inputs are two, sampling a brand new measure (only used during initialization)
% If the number of inputs are all five, sampling a new value for only one element of the measure
% This code samples a new measure value from truncated Gaussian distribution.
%
% INPUT
%   nSources - number of sources
%   lowerindex - the cell that stores all the corresponding subsets (lower index) of measure elements
%   upperindex - the cell that stores all the corresponding supersets (upper index) of measure elements
%   prevSample (optional)- previous measure, before update
%   IndexsubsetToChange (optional)- the measure element index to be updated
%   sampleVar (optional) - sampleVar set in Parameters
% OUTPUT
%   - measure - new measure after update
%
% Written by: X. Du 03/2018
%
z = rand(1);

if z < eta && nargin > 4
    
    %get upper/lower bounds
    Nmeasure = 2^nSources-1;
    measure = prevSample;
    if (IndexsubsetToChange<=nSources) && (IndexsubsetToChange>=1) %singleton
        lowerBound = 0; 
        upperBound = min(measure(upperindex{IndexsubsetToChange})); 
    elseif (IndexsubsetToChange>=(Nmeasure-nSources)) && (IndexsubsetToChange<=(Nmeasure-1)) %(nSources-1)-tuple
        lowerBound = max(measure(lowerindex{IndexsubsetToChange})); 
        upperBound = 1;
    else  %remaining elements
        lowerBound = max(measure(lowerindex{IndexsubsetToChange})); 
        upperBound = min(measure(upperindex{IndexsubsetToChange})); 
    end
    
    %update subset in measure based on bounds
    if lowerBound == 1 && upperBound == 1
        measure(IndexsubsetToChange) = 1;
    elseif lowerBound == 0 && upperBound == 0
        measure(IndexsubsetToChange) = 0;
    elseif lowerBound == 0 && upperBound == 1
        measure(IndexsubsetToChange) = randi([0 1], 1, 1);
    end
    
else
    
    coin = rand(1);
    
    if coin >= 0.5
        %bottom up approach
        measure = zeros(1,(2^nSources-1));
        measure(1:nSources) = randi([0 1],1,nSources); %sample densities
        measure(end) = 1;
        for j = (nSources+1) : (size(measure,2)-1)
            lowerBound = max(measure(lowerindex{j}));
            if lowerBound == 1
                measure(j) = 1;
            else
                measure(j) = randi([0 1], 1, 1);
            end
        end
    else
        %top down approach
        Nmeasure = 2^nSources-1;
        measure = zeros(1, Nmeasure);
        measure(Nmeasure) = 1;
        measure(Nmeasure-1:-1:Nmeasure-nSources) = randi([0 1], 1, nSources);
        
        for j = (Nmeasure-nSources-1):-1:1
            upperBound = min(measure(upperindex{j}));
            if upperBound == 1
                measure(j) = randi([0 1], 1, 1);
            else
                measure(j) = 0;
            end
        end
    end
    
end

end

