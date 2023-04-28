#' Identify duplicate samples and remove one 
#'
#' remove.duplicate.samples() receives data in list format. 
#' Keeps sample missing least data, and thereafter, by precedence in file
#'
#' @param dart_data      -- dart data list  [required]
#' @param least_missing  -- if TRUE, by least missing   
#' @param seed           -- random seed    
#' @export
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @examples
#' remove.duplicate.samples(dart_data, least_missing=TRUE, remove_fixed_loci=TRUE)
#
 
remove.duplicate.samples_MC <- function(dart_data, least_missing=TRUE, remove_fixed_loci=TRUE) {

   cat("  Looking for duplicate sample names \n")
   rownames(dart_data$gt) <- dart_data$sample_names 
   sample_names <- rownames(dart_data$gt)

   sample_counts <- as.matrix(table(as.character(sample_names)))
   duplicated_samples <- rownames(sample_counts)[which(sample_counts > 1)]

   num_duplicated_samples <- length(duplicated_samples)
  
   if (num_duplicated_samples < 1) {
       cat(" no sample names were duplicated, nothing to do \n")
       return(dart_data)
   } else {
      cat(" found ", num_duplicated_samples," samples that are duplicated \n")
      cat(" will keep one set of genotype data for each \n")

      missing  <- is.na(dart_data$gt)
      count_of_missing_by_sample  <- rowSums(missing) 

      for (n in 1:num_duplicated_samples) {
      
         sample            <- duplicated_samples[n]
         indices_of_sample <- which(sample_names == sample)
          
         ind_min_missing <- which.min(count_of_missing_by_sample[ indices_of_sample ])
         indices_to_remove <- indices_of_sample[ -ind_min_missing ]

         if (n == 1) { indices_to_remove_cumulative <- indices_to_remove } 
         if (n > 1)  { indices_to_remove_cumulative <- c(indices_to_remove_cumulative, indices_to_remove) } 
      }

      dart_data_proc <- dart_data
      dart_data_proc$gt <- dart_data$gt[ -indices_to_remove_cumulative , ] 
      dart_data_proc$sample_names <- dart_data_proc$sample_names[ -indices_to_remove_cumulative ]
      dart_data_proc$treatment    <- paste(dart_data_proc$treatment, "dupSampRm", sep="_")

      if (remove_fixed_loci) {
         cat("  Proceeding with removal of loci that become fixed \n")
         dart_data_proc <- remove.fixed.snps(dart_data_proc)
             
      } else {
         cat("  Warning: samples removed, loci might be fixed in the remaining samples. Consider removing non-polymorphic sites \n")
      }

      return(dart_data_proc)
   }

}
 
