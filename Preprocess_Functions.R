##############
# PreProcess_Functions.R 
# R code to pre-process output from Josh Schwab (JS) or Yea-Hung Chen (YHC)
#
# Laura B. Balzer, PhD MPhil
# lbalzer@umass.edu
# Lead Statistician for SEARCH
#---------------------------

#################


#---------------------------------
# preprocess: function to take in JS or YHC data and preprocess for running the Stage 2 analysis: 
#   codes community as id
#   codes intervention as A
#   adds dummy column U=1


preprocess <- function(data.input, YHC=F){
	
	# transform pairs to be numeric
	data.input$pair <- as.numeric(as.character(data.input$pair))

	if(YHC){
		
		# rename community to be id
		names(data.input)[grep('community', names(data.input))] <- 'id'
		
		# change the intervention variable to be numeric	
		data.input <- cbind(data.input, A=as.numeric(as.logical(data.input$intervention)) )
		
		# double check that using the unblinded dataset
		load("outputs-postbl-unblinded.RData")
		check <- aggregate(outputs[,'intervention'], by=list(outputs$community_number),  mean)[,2]
		if( sum(check != data.input$A)!=0 ){
			print("STOP STOP STOP STOP")
		}
		rm(outputs, check)
	
	} else{
		
		# Exclude anyone that was flagged as an SEARCH-id related error
		data.input <- subset(data.input, !(data_flag | dead_0 | move_0) )
		
		# independent unit is community_number
		names(data.input)[grep('community_number', names(data.input))] <- 'id'
	
		data.input <- cbind(data.input, A= as.numeric(as.logical(data.input$intervention) ) )

		# double check that using the unblinded dataset
		load("outputs-postbl-unblinded.RData")
		check <- aggregate(outputs[,'intervention'], by=list(outputs$community_number),  mean)[,2]
		check2 <- aggregate(data.input$A, by=list(data.input$id), mean)[,2]
		if( sum(check != check2)!=0 ){
			print("STOP STOP STOP STOP")
		}
		rm(outputs, check, check2)
}

	# if haven't added a dummy variable for unadjusted
	if( sum(grep('U', colnames(data.input)))==0){
		data.input <- cbind(U=1, data.input) 
	}
	
	print('***preprocessing done***')
	
	data.input
}

