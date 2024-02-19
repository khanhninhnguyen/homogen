# Function to check if a file exists and add its description to the README content
# Function to check if a file exists and add its description to the README content
add_file_description <- function(description) {
  return(c(description, "\n"))
}

# File descriptions
summary_list_selected <- c(
  "list_selected_nmin200_10nearby.txt", 
  "Examines change-points with up to 10 nearby stations, filtered by distance and length. Structure:",
  "- Column 1: Main - Name of the main station.",
  "- Column 2: Brp - Detected change-point.",
  "- Column 3: Nearby - Name of the nearby station.",
  "- Column 4: Coincident - Indicates similar change-point detection nearby within 10 days (1 = yes, 0 = no).",
  "- Column 5: Dist_noise - Distance to other change-points in the main station, flagged if ≤80 days.",
  "- Column 6: Noise - Presence in a cluster (0 = no, 1 = yes).",
  "- Columns 7-10: Main/Nearby_beg_new & Main/Nearby_end_new - Start and end days of the main/nearby series used for the test of change in mean.",
  "- Columns 11-16: N_main/nearby/joint_bef & N_main/nearby/joint_aft - Point counts in main, nearby, and joint series segments before and after the change-point.",
  "- Columns 17-22: R_main/nearby/joint_bef & R_main/nearby/joint_aft - Point rate to total length of main, nearby, and joint series segments before and after the change-point.",
  "- Columns 23-24: Dd & Dh - Horizontal (km) and vertical (m) distances.",
  "- Column 25: N_join_min - Minimum segment length of joint series.",
  "- Column 26: change - Result have an error (0 = no, 1 = yes)."
)

FGLS_jump_tvalue <- c(
  "FGLS_jump_tvalue.txt: Estimated jumps and their t-values for each test from the FGLS procedure."
)

FGLS_arma_coef <- c(
  "FGLS_arma_coef.txt: Estimated coefficients of the noise model for each series.",
  "Each series has 2 values: one for phi and one for theta.",
  "If the noise model is AR(1): Theta is NA.",
  "If the noise model is MA(1): Phi is NA.",
  "The G-E series may have two NA values; it is run only once for all nearby."
)

FGLS_MWvar <- c(
  "FGLS_MWvar.txt: Mean of the estimated moving window variance for each series."
)

FGLS_other_estimates <- c(
  "FGLS_other_estimates.txt: Other estimates from the FGLS procedure, including 8 Fourier components and the global mean."
)


# Initialize README content
readme_content <- c(
  "Attribution\n\n",
  "This package is used to attribute the change-points in relative homogenization.\n",
  "Output File Descriptions\n"
)

# Add file descriptions
readme_content <- c(readme_content, 
                    add_file_description(summary_list_selected),
                    add_file_description(FGLS_jump_tvalue),
                    add_file_description(FGLS_arma_coef),
                    add_file_description(FGLS_MWvar),
                    add_file_description(FGLS_other_estimates))

# Write the README content to a file
writeLines(readme_content, "README.txt")

# Print a message
message("README.txt template with output file descriptions has been created.")



