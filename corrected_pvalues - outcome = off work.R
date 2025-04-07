dir = "filepath"

### list files in working directory
file_list <- list.files(dir, "coeffs_")
path_list <- paste0(dir, "\\", file_list)
name_list <- gsub("coeffs_", "", gsub(".csv", "", file_list))

### set up empty list to store results
dat_list <- as.list(NULL)

### loop over files
for(i in 1:length(path_list)) {
  
  ### select model output file to read in
  path <- path_list[i]
  name <- name_list[i]
  
  ### read in model output
  dat <- read.csv(path)
  
  ### restrict model output to coefficients for current LC
  dat <- dat[grep("currlc", dat$X),]
  
  ### restrict model output to coefficients for interactions
  dat <- dat[grep(":", dat$X),]
  
  ### calculate corrected p-values
  dat$p_bonferroni <- p.adjust(dat[["Pr...z.."]], method="bonferroni")
  dat$p_holm <- p.adjust(dat[["Pr...z.."]], method="holm")
  dat$p_BH <- p.adjust(dat[["Pr...z.."]], method="BH")
  dat$p_BY <- p.adjust(dat[["Pr...z.."]], method="BY")
  
  ### add modifer variable
  dat$modifier <- name
  
  ### restrict dataset to variables of interest
  dat <- dat[,c("modifier", "X", "Pr...z..", "p_bonferroni", "p_holm", "p_BH", "p_BY")]
  
  ### rename variables
  colnames(dat)[2:3] <- c("interaction", "p_orig")
  
  ### reset row names
  rownames(dat) <- NULL
  
  ### store output dataset
  dat_list[[i]] <- dat
  
}

### stack output datasets
dat_out <- dat_list[[1]]
if(length(dat_list)>1) {
  for(i in 2:length(dat_list)) {dat_out <- rbind(dat_out, dat_list[[i]])}
}

### write out stacked dataset
write.csv(dat_out, file=paste0(dir, "\\corrected_pvalues.csv"), row.names=FALSE)
