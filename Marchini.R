

###############################################################################################
####   Implementation of a slightly modified method from "A new statistical approach to    ####
####   detecting significant activation in functional MRI" of Marchini and Ripley (2001).  ####
####   Written by:  Wiktor Olszowy, University of Cambridge                                ####
####   Contact:     wo222@cam.ac.uk                                                        ####
####   Created:     March-June 2016                                                        ####
###############################################################################################


#fMRI:            4D array fMRI dataset, 4th dimension should refer to time
#                 .nii file can be transformed to an array using function f.read.nifti.volume {AnalyzeFMRI}
#ROI:             1D vector, use only when fMRI=NULL
#wins_thr:        Winsorization threshold (applied after data normalisation -> mean 0, sd 1); for no Winsorization specify 'Inf'
#detr:            whether to detrend (using the running-lines smoother)
#TR:              the repetition time [sec]
#no_cores:        the number of cores, e.g. detectCores()
#fund_freq:       the experiment frequency [Hz=1/s]
#tr_x_axis_part:  left part of the x axis to be transformed to allow the spline more flexibility at lower frequencies (specified as proportion)
#tr_x_axis_dist:  factor with which the distances between the first 'tr_x_axis_part' of all frequencies are multiplied to allow the spline more flexibility
#
#
#return:          array containing the periodogram values and the p-values


Marchini = function(fMRI=NULL, ROI=NULL, wins_thr=3, detr=F, TR=2, no_cores=16, fund_freq=1/32, tr_x_axis_part=0.3, tr_x_axis_dist=1) {
	
	#the function works both for all voxel analysis, as well as for a ROI analysis
	#for a ROI an auxiliary 1*1*1*(time span) fMRI dataset constructed
	if (!is.null(ROI)) {
		fMRI = array(-1, dim=c(1,1,1,length(ROI)))
		fMRI[1,1,1,] = ROI
	}
	#the considered Fourier frequencies divided by TR
	freq = spec.pgram(fMRI[1,1,1,], plot=F)$freq/TR
	f_no = length(freq)
	#parameter 'span' = parameter 'k' in Marchini paper / length(data vector)
	time = (1:dim(fMRI)[4])*TR
	span = (2/fund_freq)/max(time)
	#extreme frequencies 0 and 1/(2*scanning time), as well as the fundamental frequency and its first
	#2 harmonics not included in the smoothing procedure
	fund_freq_at = which(abs(freq-fund_freq)==min(abs(freq-fund_freq)))[1]
	harm_1_at = which(abs(freq-2*fund_freq)==min(abs(freq-2*fund_freq)))[1]
	harm_2_at = which(abs(freq-3*fund_freq)==min(abs(freq-3*fund_freq)))[1]
	ind = (1:f_no)[-c(1, fund_freq_at, harm_1_at, harm_2_at, f_no)]
	#transforming the x axis: to allow the spline more flexibility at lower Fourier frequencies
	#multiplying the distances between the Fourier frequencies by 'tr_x_axis_dist' for the first 'tr_x_axis_part' of all the Fourier frequencies
	f_no_part = round(f_no*tr_x_axis_part)
	freq_tr = freq
	freq_tr[1:f_no_part] = freq_tr[1:f_no_part]*tr_x_axis_dist
	freq_tr[(f_no_part+1):f_no] = freq_tr[(f_no_part+1):f_no]+(freq[1])*f_no_part*(tr_x_axis_dist-1)
	x_freq = freq_tr[ind]
	#parallel computation
	spec_and_pval_res = simplify2array(mclapply(1:dim(fMRI)[1], function(i) {
		spec_and_pval = array(-1, dim=c((dim(fMRI)*c(1,1,1,0)+c(0,0,0,2))[2:4],f_no+1))
		for (j in 1:dim(fMRI)[2]) {
			for (k in 1:dim(fMRI)[3]) {
				data = fMRI[i,j,k,]
				#checking whether data does not only consist of zeros
				if (abs(mean(data))>0.0001 && abs(var(data))>0.0001) {
					#standardizing the data (to mean=0 and sd=1)
					data = scale(data)
					#Winsorizing the data
					data[which(data < (-1)*wins_thr)] = (-1)*wins_thr
					data[which(data > wins_thr)] = wins_thr
					if (detr==T) {
						#running-lines smoother
						supsmu_fit = supsmu(x=time, y=data, span=span)$y
						data = data - supsmu_fit
					}
					#calculating the periodogram
					P = spec.pgram(data, detrend=F, plot=F)
					#saving the periodogram
					spec_and_pval[j,k,1,1:f_no] = P$spec
					#the R test statistic
					spec_and_pval[j,k,2,2:(f_no-1)] = exp(predict(gam((log(P$spec[ind])-digamma(1)) ~
						s(x_freq)), newdata=data.frame(x_freq = freq_tr[2:(f_no-1)])))
					#the p-value
					spec_and_pval[j,k,1,f_no+1] = 1 - pexp(spec_and_pval[j,k,1,fund_freq_at]/(spec_and_pval[j,k,2,fund_freq_at]))
				}
			}
		}
		return(spec_and_pval)
	}, mc.cores=no_cores))
	#changing the order of the dimensions, since parallelisation causes the first dimension to be put at the end
	spec_and_pval = aperm(spec_and_pval_res, c(5,1,2,3,4))
	return(spec_and_pval)

}

