### Welcome!
***

Read the manuscript [here](https://www.biorxiv.org/content/10.1101/759498v4)

* ```datasets``` contains:
	* ```NR_LH_Protein_Martin``` A dataset of non-redundant ab-ag complexes from [AbDb](http://bioinf.org.uk/abs/abdb/) 
	* ```3did``` A dataset of protein-protein interactions from [3did](https://3did.irbbarcelona.org/)

* ```src``` contains codes (in ```Python``` and ```R```) for main and supplementary figures
* ```figures``` contains main and supplementary figures
* ```dl``` contains codes for deep learning
* ```sl``` contains codes for shallow learning
* ```suplemental_text``` contains codes and tex files the derivation of an analytical solution for the theoretical number of motifs within one CDR/FR region
* ```ramachandran_files``` contains codes and figures for Ramachandran angle calculations
* Analyses were primarily done using the files:
	* respairs_segment_notationx* (AbDB)
	* threedid_no_iglike_notationx* (3did)
* Raw deep learning models and outputfiles are also available at ```https://archive.sigma2.no``` under DOI ```10.11582/2020.00060```

* Dataset details

Dataset title: Structural interaction motifs data files.

Date: 2019 greifflab.org

### preprocessed
***

1. File name:

          respairs_segment_notationx_len_merged.csv
This file contains residues pairs annotated by region (segment). Along with the structural interaction motifs for paratopes and epitopes. 

          Attributes:
              pdbid: a unique identifier of an entry in Protein Data Bank (PDB).
              abchain: antibody chains (light [L] and heavy [H]).
              segment: antibody regions [FR1–3 and CDR1–3] with chain annotations. Follows the Martin numbering scheme.
              paratope: interacting residues in a paratope.
              plen: the number of residues in a paratope.
              shiftset: a set containing the residue number differences in a paratope. 
              gapset: a set containing the number of non-interacting residues (gaps) in a paratope.
              abresnumiset: a set containing the residue number in a paratope.
              ab_motif: structural interaction motif of a paratope.
              absegment: as in `segment`, without chain annotations.
              gapstatus: gap status (0 or 1).
              gapstrstatus: gap status in string (continuous or discontinuos).
              ab_motiflen: the length of structural interaction motif (paratope)
              ag_motiflen: the length of strctural interaction motif (epitope)
              epitope: interacting residues in an epitope.
              epitope_len: the number of residues in an epitope.
              ag_motif: structural interaction motif of an epitope.
              agresnumiset: a set containing the residue number in a paratope. 
              agchain: the chain of the antigen

2. File formats: comma separated file (CSV).
3. Versioning: All changes to this dataset may be documented in a changelog in this README document.
4. A number of files in this directory were derived from this file.


### dl/dataset
***

1. File name:

          motif*.tsv
          motif*pos*.tsv
          paraepi.tsv
          epipara.tsv
These files contain ~5000 pairs of paratope-epitope structural interaction motifs (motif*.tsv) [`pos` is with position] or pairs of paratope-epitope sequences (paraepi.tsv and epipara.tsv) 

          Attributes:
              epipara.tsv: first and second columns are epitope and paratope sequences respectively.
              paraepi.tsv: first and second columns are paratope and epitope sequences respectively.
              motif_epiparadash.tsv: first and second columns are epitope and paratope interaction motifs respectively.
              motif_paraepidash.tsv: first and second columns are paratope and epitope interaction motifs respectively.
              motif_epiparadash_pos.tsv: first and second columns are epitope and paratope interaction motifs respectively. With position annotation.
              motif_paraepidash_pos.tsv: first and second columns are paratope and epitope interaction motifs respectively. With position annotation.
                           
2. File formats: tab separated file (TSV).
3. Versioning: All changes to this dataset may be documented in a changelog in this README document.

### dl/dataset_ppi
***

1. File name:

          motif*.tsv
          motif*pos*.tsv
          paraepi.tsv
          epipara.tsv
These files contain ~20000 pairs of paratope-epitope structural interaction motifs (motif*.tsv) [`pos` is with position] or pairs of paratope-epitope sequences (paraepi.tsv and epipara.tsv). **For protein-protein interaction (PPI), a motif and its interacting partner are the analog to paratope and epitope in antibody-antigen scenario.** 

          Attributes:
              epipara.tsv: first and second columns are epitope and paratope sequences respectively.
              paraepi.tsv: first and second columns are paratope and epitope sequences respectively.
              motif_epiparadash.tsv: first and second columns are epitope and paratope interaction motifs respectively.
              motif_paraepidash.tsv: first and second columns are paratope and epitope interaction motifs respectively.
              motif_epiparadash_pos.tsv: first and second columns are epitope and paratope interaction motifs respectively. With position annotation.
              motif_paraepidash_pos.tsv: first and second columns are paratope and epitope interaction motifs respectively. With position annotation.
                           
2. File formats: tab separated file (TSV).
3. Versioning: All changes to this dataset may be documented in a changelog in this README document.


### NR\_LH\_Protein_Martin
***

1. File name:

          *.pdb
          
These files contain atomic coordinates (atoms and residues) of an antibody-antigen complex in PDB format.
[Go to PDB format specification.](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html)

          Attributes:
          
                           
2. File formats: protein data bank (PDB).
3. Versioning: All changes to this dataset may be documented in a changelog in this README document.

### 3did
***

1. File name:

          3did_flat.txt
          
The file contains residues pairs (protein-protein interaction, PPI) of all protein complex in pdb `pdb version 2019_01`.
[Go to 3did file specification.](https://3did.irbbarcelona.org/download.php#flat_files)

          Attributes:
          
                           
2. File formats: 3did flat file.
3. Versioning: All changes to this dataset may be documented in a changelog in this README document.


2019 GreiffLab


