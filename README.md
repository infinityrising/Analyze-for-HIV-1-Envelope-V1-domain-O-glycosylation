# Analyzer-of-HIV-1-Env-V1-domain-O-glycosylation

This program is designed to take a ".fasta" protein alignment of HIV-1 Envelope sequences and analyze them for O-glycosylation sites within the Envelope V1 domain. The output of the program is a ".xlsx" file that can be opened in Microsoft Excel.

To run this program, two files are needed:
1. A protein alignment of HIV-1 Envelope sequences
2. The output from the NetOGlyc4.0 prediction server (http://www.cbs.dtu.dk/services/NetOGlyc/). If working with a large protein alignment, predictions can be generated with the NetOGlyc4.0-BulkAnalyzer (https://github.com/infinityrising/NetOGlyc4.0-BulkAnalyzer)

Required steps:
1. Line 1: change to the name of your ".fasta" file
2. Line 2: change to the location of the ".fasta" file
3. Line 91: determine the boundaries of the V1 domain in your protein alignment. If it begins at amino acid 239 and ends at 310, select 239:310.
4. Line 149: change the name to your file containing the NetOGlyc4.0 data.

Output:
The output will be a table of values. These values will include:
SEQUENCE NAME || V1 LENGTH || # OF S/T IN V1 || # OF SCORED S/T IN V1 (>0.5) || MEAN SCORE (STD DEV) || AMINO ACID # || SCORE
