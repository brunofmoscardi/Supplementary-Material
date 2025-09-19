# Set working directory -> to source file location (optional in RStudio)

library(MixSIAR)

# ------------------------------------------------------
# Load mixture/consumer data (human isotopes)
# ------------------------------------------------------
mix <- load_mix_data(filename="Isotopes_NW_Nqn_origen_NEnqn.csv",
                     iso_names=c("d13C","d15N"), # isotope variables
                     factors="id",                # grouping factor
                     fac_random=FALSE,            # fixed effect
                     fac_nested=NULL,
                     cont_effects=NULL)
mix

# ------------------------------------------------------
# Load source data (potential dietary resources)
# ------------------------------------------------------
source <- load_source_data(filename="source_NW_Nqn_NEnqn.csv",
                           source_factors=NULL, 
                           conc_dep=FALSE,          # no concentration dependence
                           data_type="means",       # source data provided as means
                           mix)
source

# ------------------------------------------------------
# Load discrimination/TDF data (trophic enrichment factors)
# ------------------------------------------------------
discr <- load_discr_data(filename="discrimination_NW_Nqn_a_NEnqn.csv", mix)
discr

# ------------------------------------------------------
# Create isospace plot (mixture, sources, TDFs)
# ------------------------------------------------------
plot_data(filename="isospace_plot",
          plot_save_pdf=TRUE,
          plot_save_png=FALSE,
          mix, source, discr)

# ------------------------------------------------------
# Informative prior (based on zooarchaeological knowledge)
# ------------------------------------------------------
plot_prior(alpha.prior= c(0.05,1.49,1.49,0.05,1.91), source) 

# ------------------------------------------------------
# Write the JAGS model file (single individual structure)
# ------------------------------------------------------
model_filename <- "MixSIAR_model_NW_Nqn_prior_to_NEnqn.txt"
resid_err <- FALSE     # no residual error
process_err <- TRUE    # include process error
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# ------------------------------------------------------
# Run model (Bayesian MCMC, very long run)
# ------------------------------------------------------
jags.3 <- run_model(run="very long",
                    mix, source, discr, model_filename,
                    alpha.prior = c(0.05,1.49,1.49,0.05,1.91)) 

# ------------------------------------------------------
# Output options (diagnostics, summary statistics, plots)
# ------------------------------------------------------
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_NW_Nqn_verylong_priors_to_NEnqn",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,    # Gelman-Rubin diagnostic
                       heidel = FALSE,   # Heidelberger-Welch diagnostic
                       geweke = TRUE,    # Geweke diagnostic
                       diag_save = FALSE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       options(max.print = 99999),
                       diag_save_ggmcmc = FALSE)

# ------------------------------------------------------
# Process output: diagnostics, summary statistics, posterior plots
# ------------------------------------------------------
output_JAGS(jags.3, mix, source, output_options)


