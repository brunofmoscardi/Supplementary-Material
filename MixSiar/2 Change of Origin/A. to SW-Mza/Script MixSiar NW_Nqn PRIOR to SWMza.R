# Set working directory -> to source file location (optional in RStudio)

library(MixSIAR)

# ------------------------------------------------------
# Load mixture/consumer data (human isotopes)
# ------------------------------------------------------
mix <- load_mix_data(filename="Isotopes_NW_Nqn_NEW_origen_SWMza.csv",
                     iso_names=c("d13C","d15N"), # isotope variables
                     factors="id",                # grouping factor
                     fac_random=FALSE,            # fixed effect
                     fac_nested=NULL,
                     cont_effects=NULL)
mix

# ------------------------------------------------------
# Load source data (dietary resources for SW Mza)
# ------------------------------------------------------
source <- load_source_data(filename="source_SW_Mza.csv",
                           source_factors=NULL, 
                           conc_dep=FALSE,          # no concentration dependence
                           data_type="means",       # source data provided as means
                           mix)
source

# ------------------------------------------------------
# Load discrimination/TDF data (trophic enrichment factors)
# ------------------------------------------------------
discr <- load_discr_data(filename="discrimination_SW_Mza.csv", mix)
discr

# ------------------------------------------------------
# Create isospace plot (mixture, sources, and TDFs)
# ------------------------------------------------------
plot_data(filename="isospace_plot",
          plot_save_pdf=TRUE,
          plot_save_png=FALSE,
          mix, source, discr)

# ------------------------------------------------------
# Informative prior (based on zooarchaeological knowledge)
plot_prior(alpha.prior= c(0.28,1.55,1.63,0.03,0.41,2.16,0.94), source) 

# ------------------------------------------------------
# Write the JAGS model file
# ------------------------------------------------------
model_filename <- "MixSIAR_model_NW_Nqn_prior_SW_Mza.txt"
resid_err <- FALSE     # no residual error
process_err <- TRUE    # include process error
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# ------------------------------------------------------
# Run model (MCMC sampling)
# ------------------------------------------------------

# Very long run: increases convergence
jags.3 <- run_model(run="very long",
                    mix, source, discr, model_filename,
                    alpha.prior = c(0.28,1.55,1.63,0.03,0.41,2.16,0.94)) 

# ------------------------------------------------------
# Output options (summary statistics, diagnostics, plots)
# ------------------------------------------------------
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics_NW_Nqn_verylong_priors_SW_Mza",
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
# Process output: diagnostics, summaries, posterior plots
# ------------------------------------------------------
output_JAGS(jags.3, mix, source, output_options)
