

###############################################################################################
####   Overlay plot with p-values for several slices given output of Marchini().           ####
####   Written by:  Wiktor Olszowy, University of Cambridge                                ####
####   Contact:     wo222@cam.ac.uk                                                        ####
####   Created:     March-June 2016                                                        ####
###############################################################################################


#M:               output of the Marchini function; array consisting of periodogram values and p-values
#fMRI:            4D array fMRI dataset, 4th dimension should refer to time
#                 .nii format can be transformed to an array using function f.read.nifti.volume {AnalyzeFMRI}
#slices:          the axial slice numbers, for which the lowest p-values will be plotted on the fMR image
#pval_thr:        the multiple testing adjusted significance level
#crop:            the number of voxels on the left and right sides of the axial slices that can be deleted (where no brain)
#filename:        name of the file, without file extension, will be ".pdf"


Marchini_overlay_plot = function(M, fMRI, slices=c(2, 9, 16, 23, 30), pval_thr=1e-04, crop=9, filename="Mar_mod_overlays_correct") {

	#the 1st dimension of fMRI, without the empty space on the left and right sides
	d1 = dim(fMRI)[1]-crop*2
	#the 2nd dimension of fMRI
	d2 = dim(fMRI)[2]
	#multiplying the first dimension with 'length(slices)', so that 'length(slices)' images can be put next to each other
	fMRI_combined = array(-1, dim=c(d1*length(slices)+24,d2,1))
	pval_combined = array(-1, dim=dim(fMRI_combined))
	#combining several slices
	for (i in 1:length(slices)) {
		fMRI_combined[((i-1)*d1+1):(i*d1),,] = fMRI[(crop+1):(d1+crop),,slices[i],100]
		pval_combined[((i-1)*d1+1):(i*d1),,] = M[(crop+1):(d1+crop),,slices[i],1,81]
	}
	#making the colour legend
	ic = log(seq(0.1, 1, length.out=d2))
	ic_t = (1-(ic - min(ic))/max(ic - min(ic)))*1350 + 1
	#making the colour palette: the brighter, the lower the p-value
	sign_col = heat.colors(n=2048)[ic_t]
	#making the colour legend bar
	for (i in 23:19) {
		pval_combined[dim(pval_combined)[1]-i,,1] = (ic_t-min(ic_t))/max(ic_t-min(ic_t)) * pval_thr
	}
	#producing the pdf
	pdf(paste(filename, ".pdf", sep=""), width=5*(5*(dim(fMRI)[1]-2*crop)/(dim(fMRI)[1])+24/(dim(fMRI)[1])), height=5)
		overlay(nifti(fMRI_combined), nifti(pval_combined), z=1,
			plot.type="single", zlim.y=c(0,pval_thr), col.y=sign_col)
		mtext(text=pval_thr, at=0.95, side=1, line=3.5, col="white", cex=2.6, adj=0)
		mtext(text=0, at=0.95, side=3, line=2, col="white", cex=2.6, adj=0)
	dev.off()

}

