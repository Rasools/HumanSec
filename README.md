# HumanSec
Toolbox for generating protein-secretion genome-scale model for Human.

To test this toolbox and generate a protein-specific genome-scale model you just need to run humanSec.m function. This function needs following inputs:
  1) A GEM model (Human1 or Human1-drived models) specialized to support protein secretory reactions. you can always find  latest version of this model in input folder named as ihumanSec.
  2) Uniprot ID for your protein.
  3) Protein-Specific Information Matrix (there is a list of human proteins in excel format in "inputs" folder named PSIM_HUMAN.xlsx). You can call this list using importPsimFromExcel.m function and then use output of this function in humanSec.m function.
  4) List of template reactions (there is a list of template reactions in excel format in "inputs" folder named tempRxns_v02.xlsx). You can call this list using importTempRxnsFromExcel.m function and then use output of this function in humanSec.m function.
