# HumanSec
Toolbox for generating protein-secretion genome-scale model for Human.

To test this toolbox and generate a protein-specific genome-scale model you just need to run humanSec.m function. This function needs following inputs:
  1) A model (HMR2 or HMR2-drived models)
  2) Uniprot ID for your interested protein
  3) Protein-Specific Information Matrix (there is a list of human proteins in excel format in repository named PSIM_HUMAN.xlsx). you can call this list using importPsimFromExcel.m function and then use it in humanSec.m function.
  4) List of template reactions (there is a list of template reactions in excel format in repository named tempRxns_v01.xlsx). you can call this list using importTempRxnsFromExcel.m function and then use it in humanSec.m function.
