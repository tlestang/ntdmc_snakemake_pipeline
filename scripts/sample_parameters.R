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

    prevalence_output <- read.csv(output_file)
    ## TODO What does the wrapped model return?
    return(100 * prevalence_output[, dim(prevalence_output)[2]])
}

prevalence_map = read.csv(snakemake@input[[2]])
ess_not_reached <- FALSE
param_and_weights <- withCallingHandlers(
    trachomAMIS::amis(prevalence_map = prevalence_map,
                      transmission_model = wrapped_model,
                      snakemake@params
                      ),
    warning = function(e) ess_not_reached <<- TRUE
)
write.csv(param_and_weights, snakemake@output[[1]], row.names=FALSE)
