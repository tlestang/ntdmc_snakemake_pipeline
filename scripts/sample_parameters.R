reticulate::use_condaenv("trachoma", required = TRUE)
trachoma_module <- reticulate::import("trachoma")
model_func <- trachoma_module$Trachoma_Simulation

wrapped_model <- function(seeds, beta_values) {
    ## write input on disk
    input_file=sprintf("beta_values_job%s.csv", snakemake@wildcards[["LEVEL"]])
    output_file=sprintf("prevalence_job%s.csv", snakemake@wildcards[["LEVEL"]])

    write.csv(
        cbind(seeds, beta_values),
        file=input_file,
        row.names=FALSE
    )

    model_func(
        input_file,
        snakemake@input[[1]],
        output_file,
        sprintf("infection_job%s.csv.csv", snakemake@wildcards[["LEVEL"]]),
        SaveOutput = F,
        OutSimFilePath = NULL,
        InSimFilePath = NULL
    )

    # "prev" == "prevalence"
    prev_output <- read.csv(output_file)

    ## Trachoma model output file is
    ##            seed, beta, time1, time2, time3 ...
    ## sample 1
    ## sample 2
    ## ...
    ## We are only interested in the last column
    end_prev_values = 100 * prev_output[, dim(prev_output)[2]]
    ## Return value should be nsamples x ntimepoints (1) matrix
    returned_prev = matrix(end_prev_values, length(end_prev_values), 1)
    return(returned_prev)
}

# 1D exponential prior function
rprior <- function(n) {
  params<-matrix(NA,n,1)
  colnames(params)<-c("beta")
  params[,1]<-rexp(n)
  return(params)
}

dprior <- function(x,log=FALSE) {
  if (log) {
    return(sum(dexp(x,log=T)))
  } else {
    return(prod(dexp(x)))
  }
}
prior<-list(rprior=rprior,dprior=dprior)

prevalence_map = read.csv(snakemake@input[[2]])

param_and_weights = trachomAMIS::amis(prevalence_map = prevalence_map,
                                      transmission_model = wrapped_model,
                                      prior = prior,
                                      amis_params = snakemake@params,
                                      seed = 1
                                      )

write.csv(param_and_weights, snakemake@output[[1]], row.names=FALSE)
