# fusion-visualization

This is a tool intended to help visualize gene fusions detected from longread RNA sequencing. Currently it only works with the output files from the flair-fusion pipeline. 
* This can be run through 9-28-bam-to-fusions-pipe.py, with additional functionality from 
 * remapFusions.py - makes zoomed in view and improves accuracy of fusion point
 * matchIsoforms.py - matches isoforms to reads and provides user with more condensed view (must be used with isoform reads from flair collapse (https://github.com/BrooksLabUCSC/flair)
 * multiMapReads.py - maps locations of reads on chromosome back to fasta, providing additional information and confidence about whether multimapping is the result of a fusion or repetitive/common sequence

The default files loaded are from an amplicon fusion dataset, so will have more reads per fusion than most files.

Click on any row of the table to view that fusion

Click on the title row of the table to sort by any column

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
  * The number in the supporting reads column of the table represents the number of reads that double map when remapping
You can also load in an optional fusion fragments aligned to the reads file
* This will be used in the remapped reads view - if you have this file loaded, you can click on any read in this view (click only on the left plot)
* After selecting a read, two plots will be displayed on the right
  * The top plot shows a full length view of the read and which portions map to which genes. 
  * The lower plot shows a zoomed in view, with the mapping location names on the y-axis and base pair mapping to the reference sequence
  * You can click on any location on the full length view and move the zoomed in view to that location
  * To exit this view and return to viewing all reads for a fusion, click on the fusion in the table to the left
