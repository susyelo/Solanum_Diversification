
#######################################################
# Code to assist models where an evolving trait can
# influence dispersal ability
#######################################################

add_jts_to_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, numtrait_states)
	{
	jts_txt_matrix = matrix(data="", nrow=numtrait_states, ncol=numtrait_states)
	for (jts_i in 1:numtrait_states)
		{
		for (jts_j in 1:numtrait_states)
			{
			newtxt = paste0("jt", jts_i, jts_j)
			jts_txt_matrix[jts_i,jts_j] = newtxt
			
			if (jts_i != jts_j)
				{
				tmprow = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d",]
				row.names(tmprow) = newtxt
				tmprow$type = "fixed"
				tmprow$init = 0.0
				tmprow$min = 0.0
				tmprow$max = 1.0
				tmprow$est = 0.0
				tmprow$desc = paste0("prop. j events where trait changes state ", jts_i, "->", jts_j)
				BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = rbind(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table, tmprow)
				}
			} # END for (jts_j in 1:numtrait_states)
		} # END for (jts_i in 1:numtrait_states)
	
	BioGeoBEARS_run_object$jts_txt_matrix = jts_txt_matrix
	
	return(BioGeoBEARS_run_object)
	}


# Add a trait (and the relevant parameters etc.

#' BioGeoBEARS_run_object A BioGeoBEARS_run_object.
#' traits_fn A traits filename. The traits should be in the same 
#' kind of format that is used for geography data. I.e., a 2-state 
#' trait will be formatted like a 2-area geography dataset would be.
add_trait_to_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, trait_fn, block_allQs=TRUE, block_allQs_trait=TRUE)
	{
	defaults='
	# Now, load a trait dataset (trait "Fly"/"Non", for flightlessness)
	trait_fn = "flightlessness.txt"
	' # END defaults

	#######################################################
	# Load the trait
	#######################################################

	# Look at your geographic range data:
	trait = getranges_from_LagrangePHYLIP(lgdata_fn=trait_fn, block_allQs=block_allQs_trait)
	trait
	
	###################################
	# Error check on trait
	###################################
	tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=BioGeoBEARS_run_object$geogfn, block_allQs=block_allQs)
	species_from_geog = sort(row.names(tipranges@df))

	species_from_traits = sort(row.names(trait@df))
	matchTF1 = all(species_from_geog %in% species_from_traits)
	matchTF2 = all(species_from_traits %in% species_from_geog)
	
	if ((matchTF1 + matchTF2) != 2)
		{
		txt = "STOP ERROR in add_trait_to_BioGeoBEARS_run_object(). The species names in your traits file do not match the species names in your geography file. These must match exactly."
		
		cat("\n\n")
		cat(txt)
		cat("\n\n")
		cat("Printing species from geography file:")
		cat("\n\n")
		species_from_geog
		cat("Printing species from traits file:")
		cat("\n\n")
		species_from_traits
		cat("\n\n")
		stop(txt)
		} # END if ((matchTF1 + matchTF2) != 2)
	###################################
	# END Error check on trait
	###################################
	
	
	# Add trait to BioGeoBEARS_run_object
	BioGeoBEARS_run_object$trait = trait
	

	#######################################################
	# Modify the BioGeoBEARS_model_object@params_table
	#######################################################
	ntrait_states = ncol(trait@df)
	ntrait_states

	#######################################################
	# Add Pmat for trait (transition rates between trait states)
	# (parameters t12, t21)
	#######################################################
	trait_Pmat_txt = matrix(data="0", nrow=ntrait_states, ncol=ntrait_states)
	param_row = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[1,]
	BGB_trait_Pmat_params_table = NULL
	for (i in 1:ntrait_states)
		{
		for (j in 1:ntrait_states)
			{
			if (i==j)
				{
				trait_Pmat_txt[i,j] = "0"
				next()
				} else {
				diag_name = paste0("t", i, j)
				trait_Pmat_txt[i,j] = diag_name
				param_row$type = "fixed"
				param_row$init = 0.001
				param_row$est = 1
				param_row$min = 0.00000001
				param_row$max = 50
				param_row$desc = "trait transition rate"
				row.names(param_row) = diag_name
				} # END if (i==j)
			BGB_trait_Pmat_params_table = rbind(BGB_trait_Pmat_params_table, param_row)
			} # END for (j in 1:ntrait_states)
		} # END for (i in 1:ntrait_states)

	# Merge t12 and t21 into params table
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = rbind(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table, BGB_trait_Pmat_params_table)

	# Add the Pmat for the trait to the BioGeoBEARS run object
	BioGeoBEARS_run_object$trait_Pmat_txt = trait_Pmat_txt


	#######################################################
	# Model for dispersal multiplier (m) when in different trait states
	# (parameter m1, m2)
	#######################################################
	mnames = paste("m", 1:ntrait_states, sep="")
	# Start with dummy rows of the right size
	BGB_trait_model_params_table = BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table[1:ntrait_states,]
	rownames(BGB_trait_model_params_table) = mnames
	# Start all the multiplier parameters at 1
	BGB_trait_model_params_table$type = "fixed"
	BGB_trait_model_params_table$init = 1
	BGB_trait_model_params_table$est = 1

	for (i in 1:ntrait_states)
		{
		BGB_trait_model_params_table$desc[i] = paste0("trait-based dispersal rate multiplier when trait=", i-1)
		}

	# Merge m1 and m2 into params table
	BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table = rbind(BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table, BGB_trait_model_params_table)

	return(BioGeoBEARS_run_object)
	} # END add_trait_to_BioGeoBEARS_run_object <- function(BioGeoBEARS_run_object, trait_fn)


