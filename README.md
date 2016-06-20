
Written by:  Wiktor Olszowy, University of Cambridge
Contact:     wo222@cam.ac.uk
Created:     March-June 2016


In this small repository there are three functions implemented in R, to be used for a modification of the approach explained in
"A new statistical approach to detecting significant activation in functional MRI" of J. Marchini and B. Ripley (2001).

  1) Marchini.R
  2) Marchini_overlay_plot.R
  3) Marchini_pgram_plot.R

The user should install and load the R packages "AnalyzeFMRI", "EnvStats", "gam", "ggplot2", "oro.nifti", "parallel"

The input to the Marchini() function should be a 4D array, can be created from a .nii file using the f.read.nifti.volume()
function from the {AnalyzeFMRI} package.

The input to the two plot functions should be the output of the Marchini() function.

"Marchini_vs_FSL_comparison.R" is the script which was used to compare the modified Marchini method with FSL (produces plots and statistics).

The comparison was discussed on the poster "fMRI experiments: can frequency-domain methods perform better than standard FSL routines?" of Wiktor Olszowy, Guy B Williams and John Aston.

The codes were tested in R 3.2.2 and under FSL 5.0.8, only under Linux Ubuntu.
