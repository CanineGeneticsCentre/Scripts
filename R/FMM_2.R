kin_est <- function(Genotypes,NoMarkers)
{

	Genvar=apply(Genotypes, 2, var)
	Genmean=apply(Genotypes, 2, mean)
	Genostand=Genotypes
	for(i in 1:NoMarkers)
	{
  		Genostand[,i]=(Genostand[,i]-Genmean[i])/sqrt(Genvar[i])
	}
	Genostand=Genostand[,(Genvar!=0)]
	Kin=cov(t(Genostand))/2

}

fmm <- function()
{
	cat("\n")
	cat(">> YOU ARE NOW IN THE R ENVIRONMENT <<\n\n")

	cat("##########################################################\n")
	cat("# Run_fast Mixed Model v5                                #\n")
	cat("#                                                        #\n")
	cat("# This program uses the outputfile from plink2fastmixmod #\n")
	cat("#                                                        #\n")
	cat("# The output is a file of adjusted chisq values that can #\n")
	cat("# be plotted in Eigenstrat_and_FMM_plot (for example)    #\n")
	cat("#                                                        #\n")
	cat("##########################################################\n\n")

	cat("\n","What is the prefix for the PLINK map and ped files? ","\n\n")
	cat("(Remember to include the suffix you used. e.g. TS_PRA_a)\n\n") 
	cat("(The output files from plink2fastmixmod must be in this directory)\n\n") 
	prefix <- scan(what ='character',n = 1, quiet = TRUE, comment.char = "");  

	cat("\nPrefix for files is: ",prefix,"\n\n")

	geno_file = paste(prefix,"_fmm.txt",sep = "")
	pheno_file = paste(prefix,"_pheno.txt",sep = "")
	out_file = paste(prefix,"_fmm.out",sep = "")
	map_file =  paste(prefix,".map",sep = "")

	if(file.exists(geno_file) == FALSE){cat("File ",geno_file," does not exist\n");cat("\n\n");return()}
	if(file.exists(pheno_file) == FALSE){cat("File ",pheno_file," does not exist\n");cat("\n\n");return()}
	if(file.exists(map_file) == FALSE){cat("File ",map_file," does not exist\n");cat("\n\n");return()}


	cat ("Genotype file:  \t",geno_file,"\n")
	cat ("Phenotype file: \t",pheno_file,"\n")
	cat ("Map file:       \t",map_file,"\n\n")
	cat ("Output file:    \t",out_file,"\n\n\n")

	run_fastmixmod (geno_file,pheno_file,out_file,map_file)
}

run_fastmixmod <- function(genotypes_file,phenotypes_file,results_file,map_file)
{

	#############################################
	# Load the chr, snp_names and map positions #
	#############################################
	cat ("Loading data from map file ",map_file,"\n\n")
	map <- read.table(map_file)
	
	chromosome <- map [c(1)]
	snp_name <- map [c(2)]
	position <- map [c(4)]


	##################################
	# Create consecutive row indexes #
	##################################

	indexes <- row(chromosome)


	##############################################
	# Add headings for these single column files #
	##############################################

	names(chromosome) = c("CHR")
	names(snp_name) = c("SNP")
	names(position) = c("POS")
	


	###################################################
	# Load in the genotypes into the matrix Genotypes #
	###################################################
	cat("Loading genotypes...\n\n")
	Genotypes <- as.matrix(read.table(genotypes_file))
        dimG=dim(Genotypes)
        Genotypes=as.numeric(Genotypes)
        dim(Genotypes)=dimG
	NoMarkers <- ncol(Genotypes)

	cat ("No. of markers: ",NoMarkers,"\n")


	##################################
	# Run kinest to estimate kinship #
	# and put the results into 'Kin' #
	##################################
	cat ("Running Kin Estimation...\n\n")
	Genvar=apply(Genotypes, 2, var)
	Genmean=apply(Genotypes, 2, mean)
	Genostand=Genotypes
	for(i in 1:NoMarkers)
	{
  		Genostand[,i]=(Genostand[,i]-Genmean[i])/sqrt(Genvar[i])
	}
	Genostand=Genostand[,(Genvar!=0)]
	Kin=cov(t(Genostand))/2


	################################################
	# Make the object Explan the same as Genotypes #
	################################################
	Explan <- Genotypes


	###################################################
	# Load in the phenotypes into the vector Response #
	###################################################
	cat ("Loading phenotypes...\n\n")
	Response <- read.table(phenotypes_file)
        Response=as.numeric(Response[[1]])
        Response=as.vector(Response)


	#########################################
	# Load in the DLL libtwovarcomp.so.0.01 #
	#########################################
	cat ("Loading the DLL...\n\n")
	dyn.load("/opt/fastmixmodtest/libtwovarcomp.so.0.01")


	##############################################
	# Set up the function so it can be used in R #
	##############################################
	FastMixedModel<-function(Response, Explan, Kin, Covariates=NULL, nu_naught=10^(-5), gamma_naught=10^5)
	{
 		Response=as.matrix(Response)
	 	Explan=as.matrix(Explan)
	 	n=length(Response)
	 	if(is.null(Covariates))
	    	Covariates=matrix(1,n,1)
 
	 	Covariates=as.matrix(Covariates)

	 	ret=.Call("rint_flmm", as.double(Explan), as.double(Response), as.integer(dim(Explan)[1]), as.integer(dim(Explan)[2]), as.double(t(Covariates)), as.integer(dim(Covariates)[2]), as.double(t(Kin)), as.double(nu_naught), as.double(gamma_naught))
  		ret

	}



	####################################################
	# Now run this function, with the three parameters #
	####################################################
	cat ("Running Fast Mixed Model...\n\n")

	cat ("(Ignore these next few lines...)\n\n")

	Result=FastMixedModel(Response, Explan, Kin)

	cat(NoMarkers)


	######################################
	# Extract chisq column from FMM file #
	######################################

	fmm_chisq <- Result [c(2)]


	#####################################
	# Add extra columns to results_file #
	# (and add the column headers)      #
	#####################################

	results_new <- cbind(indexes,chromosome,snp_name,position,fmm_chisq)
	names(results_new) = c("INDEX","CHR","SNP","POS","FMM_CHI")


	######################################################
	# Now save the Results as a file results_file        #
	######################################################
	write.table(results_new,results_file,row.names=FALSE, quote=FALSE, sep = "\t")


	######################################################
	# Now save the FMM output as a file results_file_old #
	######################################################
	results_file_old = paste(results_file,"_old",sep = "")

	write.table(Result,results_file_old,row.names=FALSE)


	#################################
	# Tell the user it has finished #
	#################################

	cat ("\n\nDon't worry about all that gobbledy-gook.  It always does that...\n\n\n")

	cat ("\n\n")
	cat ("###########################################\n")
	cat ("# Fast Mixed Model calculations completed #\n")
	cat ("###########################################\n\n")

	cat ("Results written to the file: ",results_file,"\n\n")
	
	cat ("(This can be plotted with the Excel sheets FMM_plot or FMM_plot_HD)\n\n")

	cat ("(The simple three-column FMM data is written to the file: ",results_file_old,")\n\n")


}

fmm()
q(save="yes")