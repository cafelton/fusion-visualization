# fusion-visualization

The default files loaded are from an amplicon fusion dataset, so will have more reads per fusion than most files.
Click on any row of the table to view that fusion
When loading in files, make sure to load the top .tsv file first
If the screen freezes or stops responding to clicks, just refresh the page and it will reset to the default files and try again.

There are 3 possible file sets that you can load in:
* Original fusion and reads files: nameFusions.tsv, nameReads.bed
  * This will display the entire length of the fused read - reads may be flipped  to center the predicted breakpoint (red line)
  * If you don't see reads, don't worry, there may be small exons. Use the box zoom tool to zoom in, especially around the predicted breakpoint
  * The maximum number of reads displayed per fusion is 200, but the total number of supporting reads in the table is correct
* Isoform fusion and reads files: nameIsoformFusions.tsv, nameIsoformReads.bed
  * All of the above is true, but the spanning reads column has been updated to show the number of reads in the isoforms displayed
* Remapped fusions and reads files: nameRemappedFusions.tsv, nameRemapp-seq.bed
  * This will display the sequence near the breakpoint, with an updated breakpoint (red line)
  * Scroll to view more sequence. The total amount of sequence viewable will be equal to the -s argument that you ran remapFusions.py with (default=100)
  * The maximum number of viewable reads is 50 double-remapped and 25 single-mapped to each side
  * All reads shown are those that double mapped in the original mapping to the full genome (in original nameReads.bed)
  * The number in the supporting reads column of the table represents the number of reads that double map when remapping.
