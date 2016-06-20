

###############################################################################################
####   Frequency-domain approach of Marchini (after modification) is compared against FSL. ####
####   fMRI data of comatose patients doing a visual task (from the CRIC group) is used.   ####
####   Written by:  Wiktor Olszowy, University of Cambridge                                ####
####   Contact:     wo222@cam.ac.uk                                                        ####
####   Created:     March-June 2016                                                        ####
###############################################################################################


library(AnalyzeFMRI)
library(EnvStats) #for the Gumbel distribution
library(gam)
library(ggplot2)
library(oro.nifti)
library(parallel)

setwd("/home/wo222/Marchini_method")

source("Marchini.R")
source("Marchini_pgram_plot.R")
source("Marchini_overlay_plot.R")

#significance level
alpha = 1e-04
thres = qnorm(1-alpha)

#scan for which the figures will be made
scan_fig = 29

setwd("/home/wo222/Marchini_method/checkerboard_scans")

#function choosing the last 'n' characters from a string 'x'
substrRight = function(x, n) {
	substr(x, nchar(x)-n+1, nchar(x))
}

#selecting files that end with "rotated.nii")
files = sort(list.files()[which(substrRight(list.files(), 11)=="rotated.nii")])

#FSL 'pre-processing'
res = mclapply(1:length(files), function(i) {
	curr = files[i]
	curr_short = substr(curr, 1, nchar(curr)-12)
	system(paste("cp /home/wo222/Marchini_method/design_preproc.fsf",
		" /home/wo222/Marchini_method/checkerboard_scans/", "design_", curr_short, "_preproc.fsf" , sep=""))
	#replacing the name of the file in the 'design_*_preproc.fsf' file
	system(paste("sed -i 's/'", files[1], "'/'", curr, "'/g' ", "design_", curr_short, "_preproc.fsf", sep=""))
	system(paste("feat design_", curr_short, "_preproc.fsf", sep=""))
}, mc.cores=16)

#copying pre-processed scans from the FEAT directories
for (i in 1:length(files)) {
	curr = files[i]
	curr_short = substr(curr, 1, nchar(curr)-12)
	system(paste("cp -i /home/wo222/Marchini_method/checkerboard_scans/", substr(curr, 1, nchar(curr)-4), ".feat/",
		"filtered_func_data.nii.gz", " /home/wo222/Marchini_method/checkerboard_scans/", substr(curr, 1, nchar(curr)-4),
		"_preproc.nii.gz", sep=""))
	system(paste("gunzip ", substr(curr, 1, nchar(curr)-4), "_preproc.nii.gz", sep=""))
}

#selecting files that end with "rotated_preproc.nii")
files = sort(list.files()[which(substrRight(list.files(), 19)=="rotated_preproc.nii")])

#variables to save the numbers of significant voxels for different scans
Mar_correct_no = Mar_control_no = FSL_correct_no = FSL_control_no = rep(-1, length(files))

#creating arrays to save the modified Marchini function output for the correct experiment frequency and for a control one
Mar_correct = array(-1, dim=c(length(files),64,64,32,2,81))
Mar_control = array(-1, dim=c(length(files),64,64,32,2,81))

#performing the modified Marchini analysis for all scans (with 16 cores around 1h-1.5h)
system.time(for (file in files) {
	cat(file, "\n")
	fMRI = f.read.nifti.volume(paste("/home/wo222/Marchini_method/checkerboard_scans/", file, sep=""))
	Mar_correct[which(file==files),,,,,] = Marchini(fMRI=fMRI, detr=F, tr_x_axis_part=0.3, tr_x_axis_dist=1, fund_freq=1/32)
	Mar_control[which(file==files),,,,,] = Mar_correct[which(file==files),,,,,]
	#Fourier frequency 16 = 1/(10+10)
	Mar_control[which(file==files),,,,1,81] = 1-
		pexp(Mar_control[which(file==files),,,,1,16]/Mar_control[which(file==files),,,,2,16])
})

#counting significant voxels for the modified Marchini
for (scan in 1:length(files)) {
	Mar_pval_correct = Mar_correct[scan,,,,1,81]
	Mar_pval_control = Mar_control[scan,,,,1,81]
	Mar_correct_no[scan] = sum(Mar_pval_correct>=0 & Mar_pval_correct<alpha)
	Mar_control_no[scan] = sum(Mar_pval_control>=0 & Mar_pval_control<alpha)
}

#FSL 'stats' analysis
res = mclapply(1:length(files), function(i) {
	curr = files[i]
	curr_short = substr(curr, 1, nchar(curr)-12)
	system(paste("cp /home/wo222/Marchini_method/design_stats_correct.fsf",
		" /home/wo222/Marchini_method/checkerboard_scans/", "design_", curr_short, "_stats_correct.fsf" , sep=""))
	system(paste("cp /home/wo222/Marchini_method/design_stats_control.fsf",
		" /home/wo222/Marchini_method/checkerboard_scans/", "design_", curr_short, "_stats_control.fsf" , sep=""))
	#replacing the name of the file in the 'design_stats_correct.fsf' and '_control' files
	system(paste("sed -i 's/'", files[1], "'/'", curr, "'/g' ", "design_", curr_short, "_stats_correct.fsf", sep=""))
	system(paste("sed -i 's/'", files[1], "'/'", curr, "'/g' ", "design_", curr_short, "_stats_control.fsf", sep=""))
	system(paste("feat design_", curr_short, "_stats_correct.fsf", sep=""))
	system(paste("feat design_", curr_short, "_stats_control.fsf", sep=""))
}, mc.cores=16)

#counting significant voxels for FSL
FSL_correct = array(dim=c(32,64,64,32))
FSL_control = array(dim=c(32,64,64,32))
for (i in 1:length(files)) {
	curr = files[i]
	system(paste("gunzip /home/wo222/Marchini_method/checkerboard_scans/", 
		substr(curr, 1, nchar(curr)-4), ".feat/stats/zstat1.nii", sep=""))
	system(paste("gunzip /home/wo222/Marchini_method/checkerboard_scans/", 
		substr(curr, 1, nchar(curr)-4), "+.feat/stats/zstat1.nii", sep=""))
	FSL_correct[i,,,] = f.read.nifti.volume(paste("/home/wo222/Marchini_method/checkerboard_scans/", 
		substr(curr, 1, nchar(curr)-4), ".feat/stats/zstat1.nii", sep=""))[,,,1]
	FSL_control[i,,,] = f.read.nifti.volume(paste("/home/wo222/Marchini_method/checkerboard_scans/", 
		substr(curr, 1, nchar(curr)-4), "+.feat/stats/zstat1.nii", sep=""))[,,,1]
}

#calculating the number of voxels that FSL considers to be significant (for the alpha significance level)
for (scan in 1:length(files)) {
	zstat = FSL_correct[scan,,,]
	FSL_correct_no[scan] = sum(zstat >= thres)
	zstat = FSL_control[scan,,,]
	FSL_control_no[scan] = sum(zstat >= thres)
}

#common voxel proportions (modified Marchini vs. FSL)
com_vox_prop = array(-1, dim=c(length(files),2))
for (scan in 1:length(files)) {
	Mar_pval_correct = Mar_correct[scan,,,,1,81]
	Mar_pval_control = Mar_control[scan,,,,1,81]
	zstat = FSL_correct[scan,,,]
	FSL_pval_correct = 1-pnorm(zstat)
	zstat = FSL_control[scan,,,]
	FSL_pval_control = 1-pnorm(zstat)
	com_vox_prop[scan,1] = sum(Mar_pval_correct>=0 & Mar_pval_correct<=alpha & FSL_pval_correct>=0 & FSL_pval_correct<=alpha) /
		sum(FSL_pval_correct>=0 & FSL_pval_correct<=alpha)
	com_vox_prop[scan,2] = sum(Mar_pval_control>=0 & Mar_pval_control<=alpha & FSL_pval_control>=0 & FSL_pval_control<=alpha) /
		sum(FSL_pval_control>=0 & FSL_pval_control<=alpha)
}

setwd("/home/wo222/Marchini_method")

#checking whether the test statistic under H_0 of no response with frequency w_c is indeed ~Exp(1)
pdf("Marchini_Gumbel_ind.pdf")
	ind = (1:80)[-c(1:5, 10, 20, 30, 75:80)]
	d1 = data.frame(ratio = as.vector(log(Mar_correct[scan_fig,,,,1,ind]/Mar_correct[scan_fig,,,,2,ind]))[which(Mar_correct[scan_fig,,,,2,ind]>0)])
	x = seq(-6, 4, by=0.1)
	d2 = data.frame(x=x, density=devd(-x))
	print(ggplot(data=d1, aes(ratio, ..density..)) + geom_histogram(fill=I(rgb(0.9, 0.2, 0.0, 0.5)), col="black", alpha=I(.5)) +
		xlim(c(-6,4)) + geom_line(data=d2, aes(x=x, y=density), colour="black") +
		labs(x=expression(log~ R[c]), y="Relative frequency"))
dev.off()

#showing the numbers of activated voxels for all patients for the two methods for the two designs
pdf("32_scans.pdf", width=6, height=6)
	d3 = data.frame(id=as.factor(rep(1:32,4)),
		analysis=c(rep("Mar_correct_no", 32), rep("Mar_control_no", 32), rep("FSL_correct_no", 32), rep("FSL_control_no", 32)),
		number=c(Mar_correct_no, Mar_control_no, FSL_correct_no, FSL_control_no))
	print(ggplot(data=d3, aes(x=id, y=number, fill=analysis)) + geom_bar(stat="identity") +
		labs(x="Patients", y="Number of activated voxels") +
		scale_fill_manual(breaks=c("Mar_correct_no", "Mar_control_no", "FSL_correct_no", "FSL_control_no"),
			labels=c("Mod. Marchini: correct design", "Mod. Marchini: wrong design", "FSL: correct design", "FSL: wrong design"),
			name="Approach: design \n", values=c("cyan4", "cyan2", "chocolate4", "chocolate2")) +
			theme(legend.position = c(0.21, 0.84)) +
		guides(fill=guide_legend(reverse=T, ncol=1)))
dev.off()

#making a pdf with the periodogram for one of the voxels
Marchini_pgram_plot(M=Mar_correct[scan_fig,,,,,], log=T, filename="pgram_log")

#making a pdf with the p-values overlaid on fMRI for the modified Marchini method
fMRI = f.read.nifti.volume(paste("/home/wo222/Marchini_method/checkerboard_scans/", files[scan_fig], sep=""))
Marchini_overlay_plot(M=Mar_correct[scan_fig,,,,,], fMRI=fMRI, pval_thr=alpha, filename="Mar_mod_overlays_correct")

#making a pdf with the p-values overlaid on fMRI for FSL
FSL_M = Mar_correct[scan_fig,,,,,]
FSL_M[,,,1,81] = 1-pnorm(FSL_correct[scan_fig,,,])
Marchini_overlay_plot(M=FSL_M, fMRI=fMRI, pval_thr=alpha, filename="FSL_overlays_correct")

#conducting t-tests under H_0 of equal mean numbers of activated voxels
t.test(Mar_correct_no, FSL_correct_no, paired=T, alternative="two.sided")$p.value # 0.30
t.test(Mar_control_no, FSL_control_no, paired=T, alternative="two.sided")$p.value # 0.11

#means
mean(Mar_correct_no) # 334
mean(FSL_correct_no) # 165
mean(Mar_control_no) # 11
mean(FSL_control_no) # 4

#medians
median(Mar_correct_no) # 18.5
median(FSL_correct_no) # 7.5
median(Mar_control_no) # 2
median(FSL_control_no) # 0.5

