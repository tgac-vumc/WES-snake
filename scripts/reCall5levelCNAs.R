# Date: March 2016
##############################################################################################
# FUNCTION DESCRIPTION

# Re-calls (called-) QDNAseq readcounts to improve 5-level copy number distinction, Using segment values
# HZ DELETION CUTOFF : LOG2 < -0.90
# AMPLIFICATION CUTOFF : LOG2 > 1
# Also reassings probabilities of double losses, losses, gains and amplifications

# input: called QDNAseq readcounts object
# output: RE-called QDNAseq readcounts object

# more detailed:
# Segment cutoffs are based on inspection of 111 copy number profiles
	# all desired HZ deletions had segment (log2) values < -0.9
	# all desired amplifications had segment (log2) values > 1
# Probabilities are changed so that the called profiles resemble the (recalled) 5 -level CNAs.
# For quick and easy purposes all probabilities are assigned 1 or 0 and are not recalculated.

##############################################################################################
# Content:

# 1. Recall single copy losses and double copy losses (HZ deletions)
	# A. Reset the calls
	# B. Reset the probabilities

# 2. Re-call single copy gains and Amplifications:
	# A. Reset the calls
	# B. Reset the probabilities

################################################################################################
reCall5levelCNAs <- function(calledRCs){


reCalledRCs <- calledRCs

tmp 			<- reCalledRCs
tmp.calls 		<- calls(tmp)
tmp.segs 		<- log2(segmented(tmp))
tmp.probdlosses <- probdloss(tmp)
tmp.probamps 	<- probamp(tmp)
tmp.probloss 	<- probloss(tmp)
tmp.probgain 	<- probgain(tmp)

################################################
# 1. Recall single copy losses and double copy losses (HZ deletions)


# A. Reset the calls
# first make all losses (-1 and -2) -1
tmp.calls[ tmp.calls < 0 ] <- -1
# then, only assign -2 to losses with log2(segment) < -1
tmp.calls[ tmp.segs < -1 ] <- -2

# B. Reset the probabilities (all probs of called losses become 1)
# for new calls == -1:
tmp.probloss[ tmp.calls == -1 ] <- 1
tmp.probdlosses[ tmp.calls == -1 ] <- 0
# for new calls == -2:
tmp.probloss[ tmp.calls == -2 ] <- 0
tmp.probdlosses[ tmp.calls == -2 ] <- 1

################################################
# 2. Re-call single copy gains and Amplifications:

# A. Reset the calls
# first make all gains (1 and 2) 1
tmp.calls[ tmp.calls > 0 ] <- 1
# then, only assign -2 to gains with log2(segment) > 1
tmp.calls[ tmp.segs > 1 ] <- 2

# B. Reset the probabilities (all probs of called gains become 1)
tmp.probgain[ tmp.calls == 1 ] <- 1
tmp.probamps[ tmp.calls == 1 ] <- 0

tmp.probgain[ tmp.calls == 2 ] <- 0
tmp.probamps[ tmp.calls == 2 ] <- 1

################################################
# Replace the old calls with the new calls
calls(reCalledRCs) 		<- tmp.calls
probdloss(reCalledRCs)	<- tmp.probdlosses
probamp(reCalledRCs)	<- tmp.probamps
probloss(reCalledRCs)	<- tmp.probloss
probgain(reCalledRCs)	<- tmp.probgain

reCalledRCs

} # EOF
