# MIMRF-BFM
**Efficient Multi-Resolution Fusion for Remote Sensing Data with Label Uncertainty** 

_Hersh Vakharia and Xiaoxiao Du_

[[`IEEEXplore`](https://ieeexplore.ieee.org/document/10282851)]

## Installation Prerequisites
This code uses MATLAB Statistics and Machine Learning Toolbox, MATLAB Optimization Toolbox and MATLAB Parallel Computing Toolbox.  

## Demo
Run `demo_main.m` in MATLAB.

## Main Functions
The MIMRF-BFM Algorithm runs using the following function:  
`[measure, initialMeasure, Analysis] = learnCIMeasure_bfm_multires(Bags, Labels, Parameters)`

## Inputs
The _Bags_ input is a 1xNumTrainBags cell array. Inside each cell, NumPntsInBag x nSources cell. Inside each cell, the "collection" of all possible combinations generated from the multi-resolution data set.

The Labels input is a 1xNumTrainBags double vector that takes values of "1" and "0" for two-class classification problems - Training labels for each bag.

## Parameters
The parameters can be set in the following functions:  
`[Parameters] = learnCIMeasureParams()`  

The parameters is a MATLAB structure with the following fields:
1. eta: percentage of time to make small-scale mutation
2. analysis: if = "1", save all intermediate results
3. exaustiveSearchThresh: count threshold for number of repeated samples
4. fitnessUpdateThresh: count threshold for number of times new BFM samples do not improve over past iterations

## License
This source code is licensed under the license found in the [`LICENSE`](LICENSE) file in the root directory of this source tree.

This product is Copyright (c) 2023 H. Vakharia and X. Du. All rights reserved.

## Citing MIMRF-BFM
If you use the MIMRF-BFM multi-resolution fusion algorithm, please cite the following reference using the following BibTeX entries.
```
@INPROCEEDINGS{10282851,
  author={Vakharia, Hersh and Du, Xiaoxiao},
  booktitle={IGARSS 2023 - 2023 IEEE International Geoscience and Remote Sensing Symposium}, 
  title={Efficient Multi-Resolution Fusion for Remote Sensing Data with Label Uncertainty}, 
  year={2023},
  volume={},
  number={},
  pages={6326-6329},
  doi={10.1109/IGARSS52108.2023.10282851}
}
```
## Related Work

MIMRF with normal fuzzy measures: [[`arXiv`](https://arxiv.org/abs/1805.00930)] [[`Code Repo`](https://github.com/GatorSense/MIMRF)]


Multiple Instance Choquet Integral [[`arXiv`](https://arxiv.org/abs/1803.04048)] [[`Code Repo`](https://github.com/GatorSense/MICI)]
