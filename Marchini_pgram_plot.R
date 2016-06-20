

###############################################################################################
####   Plotting the periodogram together with the estimate of the spectral density given   ####
####   the results of the Marchini() function.                                             ####
####   Written by:  Wiktor Olszowy, University of Cambridge                                ####
####   Contact:     wo222@cam.ac.uk                                                        ####
####   Created:     March-June 2016                                                        ####
###############################################################################################


#M:               output of the Marchini function; array consisting of periodogram values and p values
#i,j,k:           the coordinates of a voxel, for which the periodogram should be plotted; for a ROI give (1,1,1)
#TR:              the repetition time [sec]
#fund_freq:       the experiment frequency [Hz=1/s]
#ylim:            the length of the y axis
#log:             whether the y-axis should be with a log-scale
#filename:        name of the file, without file extension, will be ".pdf"


Marchini_pgram_plot = function(M, i=32, j=8, k=16, TR=2, fund_freq=1/32, ylim=40, log=T, filename="pgram_log") {

	freq = seq(0, 1/(2*TR), length.out=(dim(M)[5]))[2:dim(M)[5]]
	f_no = length(freq)
	dd = data.frame(freq, period=M[i,j,k,1,1:f_no], smoothed=M[i,j,k,2,1:f_no])
	dd2 = data.frame(fund=c(fund_freq, f2=2*fund_freq, f3=3*fund_freq))
	dd3 = data.frame(fund = 1)
	#producing the pdf
	pdf(paste(filename, ".pdf", sep=""), width=5.3, height=3.53333)
		if (log==F) {
			#non-log-y-scale
			print(ggplot() + labs(x="Frequency [Hz]", y="Periodogram") +
				ylim(c(0,ylim)) +
				geom_vline(data=dd2, aes(xintercept=fund), lty="dashed") +
				geom_hline(data=dd3, aes(yintercept=fund, color="Fundamental frequency + first 2 harmonics"), lty="blank") +
				geom_point(data=dd, aes(x=freq, y=period, color="Periodogram")) +
				geom_line(data=dd[2:(f_no-1),], aes(x=freq, y=smoothed, color="Spectral density estimate")) +
				scale_colour_manual(
					breaks=c("Periodogram", "Spectral density estimate", "Fundamental frequency + first 2 harmonics"),
					values=c("black", "red", "blue"),
					guide = guide_legend(override.aes = list(linetype=c("blank", "solid", "dashed"), shape=c(16,NA,NA)))) +
				theme(legend.position = c(0.69, 0.84), legend.title=element_blank()))
		} else {
			#log-y-scale
			print(ggplot() + labs(x="Frequency [Hz]", y="Periodogram") +
				scale_y_log10(limits = c(min(0.0001, min(dd$period)), ylim)) +
				geom_vline(data=dd2, aes(xintercept=fund), lty="dashed") +
				geom_hline(data=dd3, aes(yintercept=fund, color="Fundamental frequency + first 2 harmonics"), lty="blank") +
				geom_point(data=dd, aes(x=freq, y=period, color="Periodogram")) +
				geom_line(data=dd[2:(f_no-1),], aes(x=freq, y=smoothed, color="Spectral density estimate")) +
				scale_colour_manual(
					breaks=c("Periodogram", "Spectral density estimate", "Fundamental frequency + first 2 harmonics"),
					values=c("black", "red", "blue"),
					guide = guide_legend(override.aes = list(linetype=c("blank", "solid", "dashed"), shape=c(16,NA,NA)))) +
				theme(legend.position = c(0.69, 0.84), legend.title=element_blank()))
		}
	dev.off()

}

