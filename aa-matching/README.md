# aa-matching
Amino Acid Matching - Module containing functions for amino acid matching



Key functionality - determining which amino acid residue is at which position

Python module can be used by many different projects / pipelines

- NN Serology - original RSNNS used `SAP.pm` - Gio later started using `aa_matching_msf.py`
- Solid Organ Transplant Outcomes - AA assignment is done in `aa_mm_biopython.py` where functions in `aa_matching.py` were initally trialed. Will want to use AA matching.
- Hematologic Disease Association - AA assignment is done by a Perl script - `AA_mapping_6loc.pl` - currently uses Perl `SAP.pm` from NMDP and only runs on NMDP servers - info on `SAP.pm` is in `immune-gene-association` repo



## IMGT/HLA Reference Data:

Multiple sequence alignment files (MSF) are provided by IMGT/HLA: 

https://github.com/ANHIG/IMGTHLA/tree/3400/msf

Our source of amino acid positions is the Gene Feature Enumeration (GFE) database developed by NMDP:

This script in gfe-db repo loads multiple sequence alignment files (MSF) from IMGT/HLA: 

https://github.com/nmdp-bioinformatics/gfe-db/blob/master/bin/get_alignments.sh

The reason we used GFE is the intron/exon boundaries are annotated here but not IMGT/HLA.

NMDP allele frequency data groups alleles identical in antigen recognition domain (ARD).

ARD exons are exon 2 and 3 of Class I (A, B, C), and exon 2 of Class II (DRB3/4/5, DRB1, DQA1, DQB1, DPA1, DPB1)



## Load BioPython SeqRecords for HLA:

`aa_matching.py` - Loads `IMGT_HLA_Full_Protein_3330.txt` that was extracted from GFEDB into a BioPython SeqRecord object.

`aa_matching_msf.py` - Loads MSF files that come directly from IMGT/HLA database to eliminate dependency on GFEDB updates - for example GFEDB does not have IMGT/HLA v3.40.0 loaded yet.

~~TODO - Have script download MSF files from IMGT/HLA Github repo with `curl`~~

	- ~~curl needs to be installed first with `brew install curl`~~

	*GB - Used the Python `requests` module instead of curl, since it seems much simpler (reference: https://stackabuse.com/download-files-with-python/)*
		- `requests` *needs to be installed first with* `pip3 install requests`



#### Problem - over time gaps are introduced in the reference allele sequence.

MSF files have gaps - GFEDB doesn't handle gaps and we can't wait for them to do it.

`IMGT_HLA_Full_Protein_3330.txt` - has no gaps in reference allele sequence HLA-A*01:01:01:01

`IMGT_HLA_Full_Protein_3370.txt` - now has gaps that need to be removed

`IMGT_HLA_Full_Protein_3390.txt` - still more gaps that need to be removed

The position indexes for amino acids are shifting over time with new IMGT/HLA database versions.

1	2	3	4	5	6
L	A	L	T	Q	T

1	2	3	3_INS1 3_INS2	4	5	6
L	A	L	     -			-		  T	Q	T





TODO - ~~Any gapped positions in reference sequence get renamed to _INS1, _INS2, etc.~~  Developing an alternative coordinate system for the inserts.

TODO - Make a dictionary containing inserted sequence for every alleles that has sequence in the gaps.

For a hypothetical allele with an insert: A*02:430. 

1	2	3	3_INS1 3_INS2	4	5	6
L	A	L	     E			W		  T	Q	T

`seq_insert_dict`  - Dictionary that stores extra inserted sequences just in case we can use them someday

`seq_insert_dict['A*02:430']['3_INS1'] = 'E'`

`seq_insert_dict['A*02:430']['3_INS2'] = 'W'`

Real example of an allele that caused an insert is A*29:126Q

When there are gaps in the reference sequence, the IMGT/HLA alignment tool adds extra columns without incrementing positions (See alignment for HLA-A)

https://www.ebi.ac.uk/cgi-bin/ipd/imgt/hla/align.cgi

```
130        140        150                      160        170 
```

Most alleles will have a "-" for the gapped positions in reference, but a few alleles will have a residue at this position. However, these inserts are usually quite rare and not well annotated (we don't know impact of inserts).

Store all sequences in a dictionary instead of an array. Can't use an array because "INS1" can't be an index - only numerical indexes are allowed.

After renaming, we can use for loops with numbers to loop through positions again.

`for (i=1; i<=6; i++)` - would loop through the six positions to load data into an array again. 

TODO - ~~Drop the INS columns before loading in SeqRecord object for now.~~
	*GB - Dropped the INS columns after loading into the SeqRecord object using a pandas dataframe*

~~TODO - Generate a `IMGT_HLA_Full_Protein_3400.txt` file with INS positions removed for every allele.~~



-------------------------


`hlaProteinOffset` - dictionary of offsets

offsets are from examining mature protein length in IMGT/HLA alignment tool

mature protein offsets differ by locus

`getMatureProteinOffset(locus)`
