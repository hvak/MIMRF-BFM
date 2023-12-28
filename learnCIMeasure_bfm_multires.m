function [measure,initialMeasure, Analysis] = learnCIMeasure_bfm_multires(Bags, Labels, Parameters,trueInitMeasure)
%% learnCIMeasure_bfm()
% This function learns a fuzzy measure for Multiple Instance Choquet Integral (MICI)
% This function works with classifier fusion
% This function uses a minmax objective function and an evolutionary algorithm to optimize.
%
% INPUT
%    InputBags        - 1xNumTrainBags cell    - inside each cell, NumPntsInBag x nSources double. Training bags data.
%    InputLabels      - 1xNumTrainBags double  -  Bag-level training labels. 
%                   If values of "1" and "0", transform into regression bags where each negative point forms a bag.
%    Parameters                 - 1x1 struct - The parameters set by learnCIMeasureParams() function.
%    trueInitMeasure (optional) - 1x(2^nSource-1) double - true initial Measureinput
%
% OUTPUT
%    measure             - 1x(2^nSource-1) double - learned measure by MICI
%    initialMeasure      - 1x(2^nSource-1) double - initial measure by MICI during random initialization
%    Analysis (optional) - 1x1 struct - records all the intermediate results if Parameters.analysis == 1 
%                 Analysis.ParentChildMeasure - the pool of both parent and child measures across iterations
%                 Analysis.ParentChildFitness - the pool of both parent and child fitness values across iterations  
%                 Analysis.ParentChildTmpIdx2 - the indices of remaining 25% child measures selected
%                 Analysis.Idxrp - Indices after randomly change between the order of the population
%                 Analysis.measurePop - measure population
%                 Analysis.fitnessPop - fitness population
%                 Analysis.measureiter - best measure for each iteration
%                 Analysis.ratioa - min(1, exp(sum(fitnessPop)-sum(fitnessPopPrev)))
%                 Analysis.ratioachild - Percentage of children measures are kepted
%                 Analysis.ratiomVal - the maximum fitness value
%                 Analysis.ratio - exp(sum(fitnessPop)-sum(fitnessPopPrev))
%                 Analysis.JumpType - small-scale mutation or large-scale mutation
%                 Analysis.ElemUpdated - the element updated in the small-scale mutation in every iteration
%                 Analysis.subsetIntervalnPop - the intervals for each measure element for each measures in the population
%
% % Written by: X. Du 03/2018
%


%%

%set up variables
Analysis = [];
nSources = size(Bags{1}, 2); %number of sources
Parameters.nSources=nSources;
eta = Parameters.eta;
% pre-compute quantities used in evalFitness
n_bags = numel(Bags);
% Parameters.nPnts = zeros(n_bags,1);
% for i=1:n_bags
%     Parameters.nPnts(i) = size(Bags{i},1);
% end
% nPntsBags = Parameters.nPnts;

[nPntsBags, ~] = cellfun(@size, Bags);

sec_start_inds = zeros(nSources,1);
nElem_prev = 0;
for j = 1:(nSources-1)
    if j == 1 %singleton        
        sec_start_inds(j) = 0;
        MeasureSection.NumEach(j) = nSources;        % set up MeasureSection (singleton)
        MeasureSection.Each{j} = [1:nSources]';
    else  %non-singleton
        nElem_prev = nElem_prev+MeasureSection.NumEach(j-1);
        sec_start_inds(j) = nElem_prev;
        MeasureSection.NumEach(j) = nchoosek(nSources,j);%compute the cumulative number of measures of each tier. E.g. singletons,2s, 3s,..
        MeasureSection.Each{j} =  nchoosek([1:nSources],j);
    end
end

MeasureSection.NumEach(nSources) = 1;%compute the cumulative number of measures of each tier. E.g. singletons,2s, 3s,..
MeasureSection.Each{nSources} =  [1:nSources];
MeasureSection.NumCumSum = cumsum(MeasureSection.NumEach);

% %Precompute differences and indices for each bag
bag_row_ids = cell(n_bags,max(nPntsBags)); %the row indices of  measure used for each bag
diffM = cell(n_bags,max(nPntsBags));
oneV = cell(n_bags,max(nPntsBags));
for i = 1:n_bags
    for j = 1:nPntsBags(i)   %-parfor
        hx = combvec(Bags{i}{j,1},Bags{i}{j,2});
        for jj = 3:nSources
            hx = combvec(hx,Bags{i}{j,jj});
        end
        hx_unique = unique(hx','rows');
        [v, indx] = sort(hx_unique, 2, 'descend');
        vz = horzcat(zeros(size(v,1), 1), v) - horzcat(v, zeros(size(v,1), 1));
        diffM{i,j} = vz(:, 2:end);
        bag_row_ids{i,j} = zeros(size(hx_unique,1),nSources-1);
        oneV{i,j} = ones(size(hx_unique,1), 1); %Create oneV cell matrix
      for ns = 2:(nSources-1) %# of sources in the combination (e.g. j=1 for g_1, j=2 for g_12)
                elem = MeasureSection.Each{ns};%the number of combinations, e.g., (1,2),(1,3),(2,3)
            for n = 1:size(hx_unique,1)
                %singleton
                bag_row_ids{i,j}(n,1) = indx(n,1);
                %Non-singleton
                temp = sort(indx(n,1:ns), 2);
                [~,~,row_id] = ismember_findrow_mex_my(temp,elem);                
                bag_row_ids{i,j}(n,ns) = sec_start_inds(ns) + row_id;
            end   
       end

    end
end




% compute lower bound index
nElem_prev = 0;
nElem_prev_prev=0;

for i = 2:nSources-1 %sample rest
    nElem = MeasureSection.NumEach(i);%the number of combinations, e.g.,3
    elem = MeasureSection.Each{i};%the number of combinations, e.g., (1,2),(1,3),(2,3)
    nElem_prev = nElem_prev+MeasureSection.NumEach(i-1);
    if i==2
        nElem_prev_prev = 0;
    elseif i>2
        nElem_prev_prev  = nElem_prev_prev + MeasureSection.NumEach(i-2);
    end
    for j = 1:nElem
        Parameters.lowerindex{nElem_prev+j}  =[];
        elemSub = nchoosek(elem(j,:), i-1);
        for k = 1:length(elemSub) %it needs a length(elemSub) because it is taking subset of the element rows
            tindx = elemSub(k,:);
            [Locb] = ismember_findRow(tindx,MeasureSection.Each{i-1});
            Parameters.lowerindex{nElem_prev+j} = horzcat(Parameters.lowerindex{nElem_prev+j} , nElem_prev_prev+Locb);
        end
    end
end

%compute upper bound index
nElem_nextsum = 2^nSources-1;%total length of measure
for i = nSources-2:-1:1 %sample rest
    nElem = MeasureSection.NumEach(i);%the number of combinations, e.g.,3
    elem = MeasureSection.Each{i};%the number of combinations, e.g., (1,2),(1,3),(2,3)
    %nElem_next = MeasureSection.NumEach(i+1);
    elem_next = MeasureSection.Each{i+1};
    
    nElem_nextsum = nElem_nextsum - MeasureSection.NumEach(i+1); %cumulative sum of how many elements in the next tier so far
    for j = nElem:-1:1
        Parameters.upperindex{nElem_nextsum-nElem+j-1}  =[];
        elemSub = elem(j,:); 
        [~,~,row_id1] = ismember_findrow_mex_my(elemSub,elem_next);
        Parameters.upperindex{nElem_nextsum-nElem+j-1} = horzcat(Parameters.upperindex{nElem_nextsum-nElem+j-1} , nElem_nextsum+row_id1-1);
    end
end



lowerindex = Parameters.lowerindex;
upperindex = Parameters.upperindex;
sampleVar = Parameters.sampleVar;

initialMeasure = sampleMeasure_bfm(eta, nSources,lowerindex,upperindex);
J_best = evalFitness_bfm_multires(Labels, initialMeasure, n_bags, nPntsBags, oneV, bag_row_ids, diffM);
measure = initialMeasure;
t = 0;
p = 0;
q = 0;
Uf = 0;

prevMeasure = measure;
P = measure;
while 1
    t = t + 1;
    if(mod(t,10) == 0)
        disp(['Iteration: ', num2str(t)]);
    end
    
    [subsetInterval] = evalInterval(prevMeasure,nSources,lowerindex,upperindex);
    sampleMultinomial_mat(subsetInterval, 1, 'descend' );
    [iinterv] = sampleMultinomial_mat(subsetInterval, 1, 'descend' );
    
    measure_t = sampleMeasure_bfm(eta, nSources, lowerindex, upperindex, prevMeasure, iinterv);
    u = 0;
    while sum(ismember(P, measure_t, 'rows')) > 0
        [subsetInterval] = evalInterval(prevMeasure,nSources,lowerindex,upperindex);
        sampleMultinomial_mat(subsetInterval, 1, 'descend' );
        [iinterv] = sampleMultinomial_mat(subsetInterval, 1, 'descend' );
        
        measure_t = sampleMeasure_bfm(eta, nSources, lowerindex, upperindex, measure_t, iinterv);
        u = u + 1;
        if u > Parameters.exaustiveSearchThresh
            Uf = 1;
            break;
        end
    end
    if Uf
        break;
    end
    p = p + 1;
    P = vertcat(P, measure_t);
    Jt = evalFitness_bfm_multires(Labels, measure_t, n_bags, nPntsBags, oneV, bag_row_ids, diffM);
    if Jt > J_best
        J_best = Jt;
        measure = measure_t;
    else
        q = q + 1;
    end
    
    if q > Parameters.fitnessUpdateThresh
        break
    end
    
    prevMeasure = measure_t;
            
end



end

