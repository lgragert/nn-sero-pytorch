## August 18th, 2020
______________________

 * Discovered offset error seemingly present in `\*AA_poly.csv` files for all loci but B (when compared to residue 
    position numbers from the nn-serology pattern files).
	* Offset seems to remain standard at two positions upstream of where it should be - maybe there were gaps in the 
	    original files?
	* Working to solve this issue by simply altering the position offsets in the `aa_matching_msf.py` script manually.
 * I wonder if it would be beneficial in the future to store this sequence information in a relational database format 
    (accessed with SQL queries). It seems like this could make it easier to add in the sequence feature annotations. On 
    the other hand, it could make the whole system overly complex. I'll do some more research on that.
 * I've switched to using PyCharm over VSCode as an IDE for now, due to the Vim emulation and feature-rich environment. 
    Not sure if I'll stay with it though.
 * Even after fixing offset issue, there are still some special cases present in regards to mismatches between original 
    polymorphisms in various HLA loci.
    * There is the distinct possibility that there were gaps present in the files that the original RSNNS generated 
        `.pat` files were based upon, which would be problematic to say the least.