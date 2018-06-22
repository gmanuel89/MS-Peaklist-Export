ms_peaklist_export <- function() {
  
  #################### MS PEAKLIST EXPORT ####################
  # All the functions in the ms_peaklist_export() function are run within the ms_peaklist_export() function's environment, but they are called from the global environment since they were previously assigned with the <<- in the functions_mass_spectrometry() function.
  # Each value that must go to the global environment is assigned with the <<- so that it can be called from any function of the ms_peaklist_export() function.
  # In the debugging phase, run the whole code block within the {}, like as if the script was directly sourced from the file.
  
  
  ### Program version (Specified by the program writer!!!!)
  R_script_version <- "2017.12.21.1"
  ### Force update (in case something goes wrong after an update, when checking for updates and reading the variable force_update, the script can automatically download the latest working version, even if the rest of the script is corrupted, because it is the first thing that reads)
  force_update <- FALSE
  ### GitHub URL where the R file is
  github_R_url <- "https://raw.githubusercontent.com/gmanuel89/MS-Peaklist-Export/master/MS%20PEAKLIST%20EXPORT.R"
  ### GitHub URL of the program's WIKI
  github_wiki_url <- "https://github.com/gmanuel89/MS-Peaklist-Export/wiki"
  ### Name of the file when downloaded
  script_file_name <- "MS PEAKLIST EXPORT"
  # Change log
  change_log <- "1. Fixed GUI\n2. Import TXT spectra allowed\n3. Peak enveloping\n4. Classwise peak filtering\n4. Added the RMS normalization\n5. Dump parameters\n6. More accurate peak enveloping\n7. Fixed focusing on Windows\n8. Allow to set tolerance in ppm\n9. Allow to generate Skyline and Median spectra"
  
  
  
  
  
  
  ############## INSTALL AND LOAD THE REQUIRED PACKAGES
  install_and_load_required_packages(c("tcltk", "parallel", "MALDIquant", "MALDIquantForeign", "XML", "foreach"), repository = NULL, update_packages = FALSE, print_messages = TRUE)
  if (Sys.info()[1] == "Windows") {
    install_and_load_required_packages("doParallel")
  } else {
    install_and_load_required_packages("doMC")
  }
  
  
  
  
  
  ###################################### Initialize the variables (default values)
  filepath_import <- NULL
  tof_mode <- "linear"
  tolerance_ppm <- 1000
  output_folder <- getwd()
  spectra_format <- "imzML"
  low_intensity_peak_removal_threshold_method <- "element-wise"
  peak_picking_algorithm <- "SuperSmoother"
  file_type_export <- "csv"
  spectra <- NULL
  peaks <- NULL
  allow_parallelization <- FALSE
  transform_data_algorithm <- NULL
  smoothing_algorithm <- "SavitzkyGolay"
  smoothing_strength <- "medium"
  baseline_subtraction_algorithm <- "SNIP"
  baseline_subtraction_algorithm_parameter <- 200
  normalization_algorithm <- "TIC"
  normalization_mass_range <- NULL
  preprocess_spectra_in_packages_of <- 0
  mass_range <- c(3000,15000)
  average_replicates <- FALSE
  averaging_method <- "mean"
  spectral_alignment_algorithm <- NULL
  spectral_alignment_reference <- "average spectrum"
  peak_deisotoping <- FALSE
  peak_enveloping <- FALSE
  peak_filtering_mode <- "whole dataset"
  parameters_matrix <- NULL
  
  
  
  
  ################## Values of the variables (for displaying and dumping purposes)
  tof_mode_value <- "Linear"
  tolerance_ppm_value <- as.character(tolerance_ppm)
  filepath_import_value <- ""
  spectra_format_value <- "imzML"
  peak_picking_algorithm_value <- "Super Smoother"
  low_intensity_peak_removal_threshold_method_value <- "element-wise"
  allow_parallelization_value <- "NO"
  transform_data_value <- "None"
  smoothing_value <- paste0("YES", "\n( ", "Savitzky-Golay", ",\n" , "medium", " )")
  baseline_subtraction_value <- paste0("YES", "\n( ", "SNIP", ",\nIterations: ", "200", " )")
  normalization_value <- paste0("YES", "\n( ", "TIC", " )")
  average_replicates_value <- "NO"
  spectral_alignment_value <- "NO"
  spectral_alignment_algorithm_value <- ""
  spectral_alignment_reference_value <- ""
  peak_deisotoping_enveloping_value <- "None"
  peak_filtering_mode_value <- "whole dataset"
  signals_avg_and_sd_value <- "Not calculated"
  mass_range_value <- "3000, 15000"
  
  
  
  
  
  
  ##################################################### DEFINE WHAT THE BUTTONS DO
  
  ##### Check for updates (from my GitHub page) (it just updates the label telling the user if there are updates) (it updates the check for updates value that is called by the label). The function will read also if an update should be forced.
  check_for_updates_function <- function() {
    ### Initialize the version number
    online_version_number <- NULL
    ### Initialize the force update
    online_force_update <- FALSE
    ### Initialize the variable that says if there are updates
    update_available <- FALSE
    ### Initialize the change log
    online_change_log <- "Bug fixes"
    # Check if there is internet connection by pinging a website
    there_is_internet <- check_internet_connection(method = "getURL", website_to_ping = "www.google.it")
    # Check for updates only in case of working internet connection
    if (there_is_internet == TRUE) {
      try({
        ### Read the file from the web (first 10 lines)
        online_file <- readLines(con = github_R_url)
        ### Retrieve the version number
        for (l in online_file) {
          if (length(grep("R_script_version <-", l, fixed = TRUE)) > 0) {
            # Isolate the "variable" value
            online_version_number <- unlist(strsplit(l, "R_script_version <- ", fixed = TRUE))[2]
            # Remove the quotes
            online_version_number <- unlist(strsplit(online_version_number, "\""))[2]
            break
          }
        }
        ### Retrieve the force update
        for (l in online_file) {
          if (length(grep("force_update <-", l, fixed = TRUE)) > 0) {
            # Isolate the "variable" value
            online_force_update <- as.logical(unlist(strsplit(l, "force_update <- ", fixed = TRUE))[2])
            break
          }
          if (is.null(online_force_update)) {
            online_force_update <- FALSE
          }
        }
        ### Retrieve the change log
        for (l in online_file) {
          if (length(grep("change_log <-", l, fixed = TRUE)) > 0) {
            # Isolate the "variable" value
            online_change_log <- unlist(strsplit(l, "change_log <- ", fixed = TRUE))[2]
            # Remove the quotes
            online_change_log_split <- unlist(strsplit(online_change_log, "\""))[2]
            # Split at the \n
            online_change_log_split <- unlist(strsplit(online_change_log_split, "\\\\n"))
            # Put it back to the character
            online_change_log <- ""
            for (o in online_change_log_split) {
              online_change_log <- paste(online_change_log, o, sep = "\n")
            }
            break
          }
        }
        ### Split the version number in YYYY.MM.DD
        online_version_YYYYMMDDVV <- unlist(strsplit(online_version_number, ".", fixed = TRUE))
        ### Compare with the local version
        local_version_YYYYMMDDVV = unlist(strsplit(R_script_version, ".", fixed = TRUE))
        ### Check the versions (from the Year to the Day)
        # Check the year
        if (as.numeric(local_version_YYYYMMDDVV[1]) < as.numeric(online_version_YYYYMMDDVV[1])) {
          update_available <- TRUE
        }
        # If the year is the same (update is FALSE), check the month
        if (update_available == FALSE) {
          if ((as.numeric(local_version_YYYYMMDDVV[1]) == as.numeric(online_version_YYYYMMDDVV[1])) && (as.numeric(local_version_YYYYMMDDVV[2]) < as.numeric(online_version_YYYYMMDDVV[2]))) {
            update_available <- TRUE
          }
        }
        # If the month and the year are the same (update is FALSE), check the day
        if (update_available == FALSE) {
          if ((as.numeric(local_version_YYYYMMDDVV[1]) == as.numeric(online_version_YYYYMMDDVV[1])) && (as.numeric(local_version_YYYYMMDDVV[2]) == as.numeric(online_version_YYYYMMDDVV[2])) && (as.numeric(local_version_YYYYMMDDVV[3]) < as.numeric(online_version_YYYYMMDDVV[3]))) {
            update_available <- TRUE
          }
        }
        # If the day and the month and the year are the same (update is FALSE), check the daily version
        if (update_available == FALSE) {
          if ((as.numeric(local_version_YYYYMMDDVV[1]) == as.numeric(online_version_YYYYMMDDVV[1])) && (as.numeric(local_version_YYYYMMDDVV[2]) == as.numeric(online_version_YYYYMMDDVV[2])) && (as.numeric(local_version_YYYYMMDDVV[3]) == as.numeric(online_version_YYYYMMDDVV[3])) && (as.numeric(local_version_YYYYMMDDVV[4]) < as.numeric(online_version_YYYYMMDDVV[4]))) {
            update_available <- TRUE
          }
        }
        ### Return messages
        if (is.null(online_version_number)) {
          # The version number could not be ckecked due to internet problems
          # Update the label
          check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdates not checked:\nconnection problems", sep = "")
        } else {
          if (update_available == TRUE) {
            # Update the label
            check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdate available:\n", online_version_number, sep = "")
          } else {
            # Update the label
            check_for_updates_value <- paste("Version: ", R_script_version, "\nNo updates available", sep = "")
          }
        }
      }, silent = TRUE)
    }
    ### Something went wrong: library not installed, retrieving failed, errors in parsing the version number
    if (is.null(online_version_number)) {
      # Update the label
      check_for_updates_value <- paste("Version: ", R_script_version, "\nUpdates not checked:\nconnection problems", sep = "")
    }
    # Escape the function
    update_available <<- update_available
    online_change_log <<- online_change_log
    check_for_updates_value <<- check_for_updates_value
    online_version_number <<- online_version_number
    online_force_update <<- online_force_update
  }
  
  ##### Download the updated file (from my GitHub page)
  download_updates_function <- function() {
    # Download updates only if there are updates available
    if (update_available == TRUE || online_force_update == TRUE) {
      # Changelog
      tkmessageBox(title = "Changelog", message = paste0("The updated script contains the following changes:\n", online_change_log), icon = "info")
      # Initialize the variable which says if the file has been downloaded successfully
      file_downloaded <- FALSE
      # Choose where to save the updated script
      tkmessageBox(title = "Download folder", message = "Select where to save the updated script file", icon = "info")
      download_folder <- tclvalue(tkchooseDirectory())
      # Download the file only if a download folder is specified, otherwise don't
      if (download_folder != "") {
        # Go to the working directory
        setwd(download_folder)
        tkmessageBox(message = paste0("The updated script file will be downloaded in:\n\n", download_folder))
        # Download the file
        try({
          download.file(url = github_R_url, destfile = paste0(script_file_name, ".R"), method = "auto")
          file_downloaded <- TRUE
        }, silent = TRUE)
        if (file_downloaded == TRUE) {
          tkmessageBox(title = "Updated file downloaded!", message = paste0("The updated script, named:\n\n", paste0(script_file_name, ".R"), "\n\nhas been downloaded to:\n\n", download_folder, "\n\nThe current window will now close and the new updated script will be loaded!"), icon = "info")
          # Destroy the window
          try(tkdestroy(window), silent = TRUE)
          # Relaunch the script
          try(source(paste0(script_file_name, ".R")), silent = TRUE)
        } else {
          tkmessageBox(title = "Connection problem", message = paste("The updated script file could not be downloaded due to internet connection problems!\n\nManually download the updated script file at:\n\n", github_R_url, sep = ""), icon = "warning")
        }
      } else {
        # No download folder specified!
        tkmessageBox(message = "The updated script file will not be downloaded!")
      }
    } else {
      tkmessageBox(title = "No update available", message = "NO UPDATES AVAILABLE!\n\nThe latest version is running!", icon = "info")
    }
    # Raise the focus on the main window (if there is)
    try(tkraise(window), silent = TRUE)
  }
  
  ### Downloading forced updates
  check_for_updates_function()
  if (online_force_update == TRUE) {
    download_updates_function()
  }
  
  ### Force check for updates
  force_check_for_updates_function <- function() {
    # Check for updates
    check_for_updates_function()
    # Display a message
    if (update_available == TRUE) {
      # Message
      tkmessageBox(title = "Update available", message = paste0("Update available!\n", online_version_number, "\n\nPress the 'DOWNLOAD UPDATE...' button to retrieve the updated script!"), icon = "info")
    } else {
      # Message
      tkmessageBox(title = "No update available", message = "No update available!", icon = "info")
    }
  }
  
  ##### Preprocessing window
  preprocessing_window_function <- function() {
    ##### Functions
    # Transform the data
    transform_data_choice <- function() {
      # Ask for the algorithm
      transform_data_algorithm_input <- select.list(c("Square root", "Natural logarithm", "Decimal logarithm", "Binary Logarithm", "None"), title = "Data transformation", multiple = FALSE, preselect = "None")
      # Raise the focus on the preproc window
      tkraise(window)
      tkraise(preproc_window)
      # Default and fix
      if (transform_data_algorithm_input == "Square root") {
        transform_data_algorithm <- "sqrt"
      } else if (transform_data_algorithm_input == "Natural logarithm") {
        transform_data_algorithm <- "log"
      } else if (transform_data_algorithm_input == "Binary Logarithm") {
        transform_data_algorithm <- "log2"
      } else if (transform_data_algorithm_input == "Decimal logarithm") {
        transform_data_algorithm <- "log10"
      } else if (transform_data_algorithm_input == "" || transform_data_algorithm_input == "None") {
        transform_data_algorithm <- NULL
      }
      # Set the value of the displaying label
      if (!is.null(transform_data_algorithm)) {
        transform_data_value <- paste0("YES", "\n( ", transform_data_algorithm_input, " )")
      } else {
        transform_data_value <- "None"
      }
      transform_data_value_label <- tklabel(preproc_window, text = transform_data_value, font = label_font, bg = "white", width = 20, height = 2)
      tkgrid(transform_data_value_label, row = 3, column = 2, padx = c(5, 5), pady = c(5, 5))
      # Escape the function
      transform_data_algorithm <<- transform_data_algorithm
      transform_data_value <<- transform_data_value
    }
    # Smoothing
    smoothing_choice <- function() {
      # Ask for the algorithm
      smoothing_algorithm_input <- select.list(c("Savitzky-Golay","Moving Average", "None"), title = "Smoothing algorithm", multiple = FALSE, preselect = "SavitzkyGolay")
      # Raise the focus on the preproc window
      tkraise(window)
      tkraise(preproc_window)
      # Default and fix
      if (smoothing_algorithm_input == "" || smoothing_algorithm_input == "Savitzky-Golay") {
        smoothing_algorithm <- "SavitzkyGolay"
      } else if (smoothing_algorithm_input == "Moving Average") {
        smoothing_algorithm <- "MovingAverage"
      } else if (smoothing_algorithm_input == "None") {
        smoothing_algorithm <- NULL
      }
      # Strength
      if (!is.null(smoothing_algorithm)) {
        smoothing_strength <- select.list(c("small", "medium", "strong", "stronger"), title = "Smoothing strength", multiple = FALSE, preselect = "medium")
        # Raise the focus on the preproc window
        tkraise(window)
        tkraise(preproc_window)
        if (smoothing_strength == "") {
          smoothing_strength <- "medium"
        }
      }
      # Set the value of the displaying label
      if (!is.null(smoothing_algorithm)) {
        smoothing_value <- paste0("YES", "\n( ", smoothing_algorithm_input, ",\n" , smoothing_strength, " )")
      } else {
        smoothing_value <- "None"
      }
      smoothing_value_label <- tklabel(preproc_window, text = smoothing_value, font = label_font, bg = "white", width = 20, height = 3)
      tkgrid(smoothing_value_label, row = 4, column = 2, padx = c(5, 5), pady = c(5, 5))
      # Escape the function
      smoothing_strength <<- smoothing_strength
      smoothing_algorithm <<- smoothing_algorithm
      smoothing_value <<- smoothing_value
    }
    # Baseline subtraction
    baseline_subtraction_choice <- function() {
      # Ask for the algorithm
      baseline_subtraction_algorithm <- select.list(c("SNIP", "TopHat", "ConvexHull", "median", "None"), title = "Baseline subtraction algorithm", multiple = FALSE, preselect = "SNIP")
      # Raise the focus on the preproc window
      tkraise(window)
      tkraise(preproc_window)
      # Default
      if (baseline_subtraction_algorithm == "") {
        baseline_subtraction_algorithm <- "SNIP"
        baseline_subtraction_algorithm_parameter <- 200
      }
      if (baseline_subtraction_algorithm == "None") {
        baseline_subtraction_algorithm <- NULL
      }
      # SNIP
      if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "SNIP") {
        baseline_subtraction_algorithm_parameter <- tclvalue(baseline_subtraction_algorithm_parameter2)
        baseline_subtraction_algorithm_parameter_value <- as.character(baseline_subtraction_algorithm_parameter)
        baseline_subtraction_algorithm_parameter <- as.integer(baseline_subtraction_algorithm_parameter)
      } else if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "TopHat") {
        baseline_subtraction_algorithm_parameter <- tclvalue(baseline_subtraction_algorithm_parameter2)
        baseline_subtraction_algorithm_parameter_value <- as.character(baseline_subtraction_algorithm_parameter)
        baseline_subtraction_algorithm_parameter <- as.integer(baseline_subtraction_algorithm_parameter)
      } else if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "median") {
        baseline_subtraction_algorithm_parameter <- tclvalue(baseline_subtraction_algorithm_parameter2)
        baseline_subtraction_algorithm_parameter_value <- as.character(baseline_subtraction_algorithm_parameter)
        baseline_subtraction_algorithm_parameter <- as.integer(baseline_subtraction_algorithm_parameter)
      }
      # Set the value of the displaying label
      if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "SNIP") {
        baseline_subtraction_value <- paste0("YES", "\n( ", baseline_subtraction_algorithm, ",\nIterations: ", baseline_subtraction_algorithm_parameter, " )")
      } else if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "TopHat") {
        baseline_subtraction_value <- paste0("YES", "\n( ", baseline_subtraction_algorithm, ",\nHalf window size: ", baseline_subtraction_algorithm_parameter, " )")
      } else if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "median") {
        baseline_subtraction_value <- paste0("YES", "\n( ", baseline_subtraction_algorithm, ",\nHalf window size: ", baseline_subtraction_algorithm_parameter, " )")
      } else if (!is.null(baseline_subtraction_algorithm) && baseline_subtraction_algorithm == "ConvexHull") {
        baseline_subtraction_value <- paste0("YES", "\n( ", baseline_subtraction_algorithm, ")")
      } else {
        baseline_subtraction_value <- "None"
      }
      baseline_subtraction_value_label <- tklabel(preproc_window, text = baseline_subtraction_value, font = label_font, bg = "white", width = 20, height = 3)
      tkgrid(baseline_subtraction_value_label, row = 5, column = 3, padx = c(5, 5), pady = c(5, 5))
      # Escape the function
      baseline_subtraction_algorithm_parameter <<- baseline_subtraction_algorithm_parameter
      baseline_subtraction_algorithm <<- baseline_subtraction_algorithm
      baseline_subtraction_value <<- baseline_subtraction_value
    }
    # Normalization
    normalization_choice <- function() {
      # Ask for the algorithm
      normalization_algorithm <- select.list(c("TIC", "RMS", "PQN", "median", "None"), title = "Normalization algorithm", multiple = FALSE, preselect = "TIC")
      # Raise the focus on the preproc window
      tkraise(window)
      tkraise(preproc_window)
      if (normalization_algorithm == "") {
        normalization_algorithm <- "TIC"
      }
      if (normalization_algorithm == "None") {
        normalization_algorithm <- NULL
      }
      # TIC
      if (!is.null(normalization_algorithm) && normalization_algorithm == "TIC") {
        normalization_mass_range <- tclvalue(normalization_mass_range2)
        normalization_mass_range_value <- as.character(normalization_mass_range)
        if (normalization_mass_range != 0 && normalization_mass_range != "") {
          normalization_mass_range <- unlist(strsplit(normalization_mass_range, ","))
          normalization_mass_range <- as.numeric(normalization_mass_range)
        } else if (normalization_mass_range == 0 || normalization_mass_range == "") {
          normalization_mass_range <- NULL
        }
      }
      # Set the value of the displaying label
      if (!is.null(normalization_algorithm) && normalization_algorithm != "TIC") {
        normalization_value <- paste0("YES", "\n( ", normalization_algorithm, " )\n")
      } else if (!is.null(normalization_algorithm) && normalization_algorithm == "TIC") {
        if (!is.null(normalization_mass_range)) {
          normalization_value <- paste0("YES", "\n( ", normalization_algorithm, ",\nrange:\n", normalization_mass_range_value, " )")
        } else {
          normalization_value <- paste0("YES", "\n( ", normalization_algorithm, " )")
        }
      } else {
        normalization_value <- "None"
      }
      normalization_value_label <- tklabel(preproc_window, text = normalization_value, font = label_font, bg = "white", width = 20, height = 4)
      tkgrid(normalization_value_label, row = 7, column = 3, padx = c(5, 5), pady = c(5, 5))
      # Escape the function
      normalization_mass_range <<- normalization_mass_range
      normalization_algorithm <<- normalization_algorithm
      normalization_value <<- normalization_value
    }
    # Spectral alignment
    spectral_alignment_choice <- function() {
      # Ask for the algorithm
      spectral_alignment_algorithm <- select.list(c("cubic", "quadratic", "linear", "lowess", "None"), title = "Spectral alignment algorithm", multiple = FALSE, preselect = "cubic")
      # Raise the focus on the preproc window
      tkraise(window)
      tkraise(preproc_window)
      # Default
      if (spectral_alignment_algorithm == "") {
        spectral_alignment_algorithm <- "None"
      }
      if (spectral_alignment_algorithm == "None") {
        spectral_alignment_algorithm <- NULL
      }
      ## Ask for the reference peaklist
      if (!is.null(spectral_alignment_algorithm)) {
        spectral_alignment_reference <- select.list(c("auto","average spectrum", "skyline spectrum"), title = "Spectral alignment reference", multiple = FALSE, preselect = "average spectrum")
        # Raise the focus on the preproc window
        tkraise(window)
        tkraise(preproc_window)
        if (spectral_alignment_reference == "") {
          spectral_alignment_reference <- "average spectrum"
        }
      } else {
        spectral_alignment_reference <- NULL
      }
      # Set the value of the displaying label
      if (!is.null(spectral_alignment_algorithm)) {
        spectral_alignment_value <- paste0("YES", "\n( ", spectral_alignment_algorithm, ",\n", spectral_alignment_reference, " )")
      } else {
        spectral_alignment_value <- "None"
      }
      spectral_alignment_value_label <- tklabel(preproc_window, text = spectral_alignment_value, font = label_font, bg = "white", width = 20, height = 3)
      tkgrid(spectral_alignment_value_label, row = 8, column = 2, padx = c(5, 5), pady = c(5, 5))
      # Escape the function
      spectral_alignment_algorithm <<- spectral_alignment_algorithm
      spectral_alignment_reference <<- spectral_alignment_reference
      spectral_alignment_value <<- spectral_alignment_value
    }
    # TOF mode
    tof_mode_choice <- function() {
      # Catch the value from the menu
      tof_mode <- select.list(c("Linear", "Reflectron"), title = "TOF mode")
      # Raise the focus on the preproc window
      tkraise(window)
      tkraise(preproc_window)
      # Default
      if (tof_mode == "" || tof_mode == "Linear") {
        tof_mode <- "linear"
      }
      if (tof_mode == "Reflectron") {
        tof_mode <- "reflectron"
      }
      # Set the value of the displaying label
      if (tof_mode == "linear") {
        tof_mode_value <- "Linear"
      } else if (tof_mode == "reflectron") {
        tof_mode_value <- "Reflectron"
      }
      tof_mode_value_label <- tklabel(preproc_window, text = tof_mode_value, font = label_font, bg = "white", width = 20)
      tkgrid(tof_mode_value_label, row = 2, column = 3, padx = c(5, 5), pady = c(5, 5))
      # Escape the function
      tof_mode <<- tof_mode
      tof_mode_value <<- tof_mode_value
    }
    # Commit preprocessing
    commit_preprocessing_function <- function() {
      # Get the values (they are filled with the default anyway)
      # Mass range
      mass_range <- tclvalue(mass_range2)
      mass_range <- as.numeric(unlist(strsplit(mass_range, ",")))
      mass_range_value <- as.character(paste(mass_range[1], ",", mass_range[2]))
      # Preprocessing
      preprocess_spectra_in_packages_of <- tclvalue(preprocess_spectra_in_packages_of2)
      preprocess_spectra_in_packages_of <- as.integer(preprocess_spectra_in_packages_of)
      preprocess_spectra_in_packages_of_value <- as.character(preprocess_spectra_in_packages_of)
      # Preprocessing
      tolerance_ppm <- tclvalue(tolerance_ppm2)
      if (tolerance_ppm == "") {
        tolerance_ppm <- NULL
        tolerance_ppm_value <- ""
      } else {
        tolerance_ppm <- as.numeric(tolerance_ppm)
        tolerance_ppm_value <- as.character(tolerance_ppm)
      }
      # Escape the function
      mass_range <<- mass_range
      mass_range_value <<- mass_range_value
      preprocess_spectra_in_packages_of <<- preprocess_spectra_in_packages_of
      preprocess_spectra_in_packages_of_value <<- preprocess_spectra_in_packages_of_value
      tolerance_ppm <<- tolerance_ppm
      tolerance_ppm_value <<- tolerance_ppm_value
      preprocessing_parameters <<- list(mass_range = mass_range, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_algorithm_parameter = baseline_subtraction_algorithm_parameter, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of, spectral_alignment_algorithm = spectral_alignment_algorithm, spectral_alignment_reference = spectral_alignment_reference)
      # Destroy the window upon committing
      tkdestroy(preproc_window)
      # Raise the focus on the preproc window
      tkraise(window)
    }
    ##### List of variables, whose values are taken from the entries in the GUI (create new variables for the sub window, that will replace the ones in the global environment, only if the default are changed)
    mass_range2 <- tclVar("")
    preprocess_spectra_in_packages_of2 <- tclVar("")
    tolerance_ppm2 <- tclVar("")
    baseline_subtraction_algorithm_parameter2 <- tclVar("")
    normalization_mass_range2 <- tclVar("")
    ##### Window
    preproc_window <- tktoplevel(bg = "white")
    tkwm.resizable(preproc_window, FALSE, FALSE)
    tktitle(preproc_window) <- "Spectral preprocessing parameters"
    #tkpack.propagate(preproc_window, FALSE)
    # Mass range
    mass_range_label <- tklabel(preproc_window, text = "Mass range", font = label_font, bg = "white", width = 20)
    mass_range_entry <- tkentry(preproc_window, textvariable = mass_range2, font = entry_font, bg = "white", width = 20, justify = "center")
    tkinsert(mass_range_entry, "end", as.character(paste(mass_range[1],",",mass_range[2])))
    # Preprocessing (in packages of)
    preprocess_spectra_in_packages_of_label <- tklabel(preproc_window, text="Preprocess spectra\nin packages of", font = label_font, bg = "white", width = 20)
    preprocess_spectra_in_packages_of_entry <- tkentry(preproc_window, textvariable = preprocess_spectra_in_packages_of2, font = entry_font, bg = "white", width = 10, justify = "center")
    tkinsert(preprocess_spectra_in_packages_of_entry, "end", as.character(preprocess_spectra_in_packages_of))
    # Tof mode
    tof_mode_label <- tklabel(preproc_window, text="Select the TOF mode", font = label_font, bg = "white", width = 20)
    tof_mode_entry <- tkbutton(preproc_window, text="Choose the TOF mode", command = tof_mode_choice, font = button_font, bg = "white", width = 20)
    # Tolerance in ppm
    tolerance_ppm_label <- tklabel(preproc_window, text="Tolerance (in ppm)", font = label_font, bg = "white", width = 20)
    tolerance_ppm_entry <- tkentry(preproc_window, textvariable = tolerance_ppm2, font = entry_font, bg = "white", width = 10, justify = "center")
    tkinsert(tolerance_ppm_entry, "end", as.character(tolerance_ppm))
    # Transform the data
    transform_data_button <- tkbutton(preproc_window, text="Data transformation", command = transform_data_choice, font = button_font, bg = "white", width = 20)
    # Smoothing
    smoothing_button <- tkbutton(preproc_window, text="Smoothing", command = smoothing_choice, font = button_font, bg = "white", width = 20)
    # Baseline subtraction
    baseline_subtraction_button <- tkbutton(preproc_window, text="Baseline subtraction", command = baseline_subtraction_choice, font = button_font, bg = "white", width = 20)
    baseline_subtraction_algorithm_parameter_entry <- tkentry(preproc_window, textvariable = baseline_subtraction_algorithm_parameter2, font = entry_font, bg = "white", width = 10, justify = "center")
    tkinsert(baseline_subtraction_algorithm_parameter_entry, "end", as.character(baseline_subtraction_algorithm_parameter))
    # Normalization
    normalization_button <- tkbutton(preproc_window, text="Normalization", command = normalization_choice, font = button_font, bg = "white", width = 20)
    normalization_mass_range_entry <- tkentry(preproc_window, textvariable = normalization_mass_range2, font = entry_font, bg = "white", width = 20, justify = "center")
    tkinsert(normalization_mass_range_entry, "end", as.character(normalization_mass_range))
    # Spectral alignment
    spectral_alignment_button <- tkbutton(preproc_window, text="Align spectra", command = spectral_alignment_choice, font = button_font, bg = "white", width = 20)
    # Commit preprocessing
    commit_preprocessing_button <- tkbutton(preproc_window, text="Commit preprocessing", command = commit_preprocessing_function, font = button_font, bg = "white", width = 20)
    ##### Displaying labels
    tof_mode_value_label <- tklabel(preproc_window, text = tof_mode_value, font = label_font, bg = "white", width = 20)
    transform_data_value_label <- tklabel(preproc_window, text = transform_data_value, font = label_font, bg = "white", width = 20, height = 2)
    smoothing_value_label <- tklabel(preproc_window, text = smoothing_value, font = label_font, bg = "white", width = 20, height = 3)
    baseline_subtraction_value_label <- tklabel(preproc_window, text = baseline_subtraction_value, font = label_font, bg = "white", width = 20, height = 3)
    normalization_value_label <- tklabel(preproc_window, text = normalization_value, font = label_font, bg = "white", width = 20, height = 4)
    spectral_alignment_value_label <- tklabel(preproc_window, text = spectral_alignment_value, font = label_font, bg = "white", width = 20, height = 3)
    #### Geometry manager
    tkgrid(mass_range_label, row = 1, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(mass_range_entry, row = 1, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(tof_mode_label, row = 2, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(tof_mode_entry, row = 2, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(tof_mode_value_label, row = 2, column = 3, padx = c(5, 5), pady = c(5, 5))
    tkgrid(transform_data_button, row = 3, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(transform_data_value_label, row = 3, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(smoothing_button, row = 4, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(smoothing_value_label, row = 4, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(baseline_subtraction_button, row = 5, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(baseline_subtraction_algorithm_parameter_entry, row = 5, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(baseline_subtraction_value_label, row = 5, column = 3, padx = c(5, 5), pady = c(5, 5))
    tkgrid(normalization_button, row = 7, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(normalization_mass_range_entry, row = 7, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(normalization_value_label, row = 7, column = 3, padx = c(5, 5), pady = c(5, 5))
    tkgrid(spectral_alignment_button, row = 8, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(spectral_alignment_value_label, row = 8, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(preprocess_spectra_in_packages_of_label, row = 9, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(preprocess_spectra_in_packages_of_entry, row = 9, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(tolerance_ppm_label, row = 10, column = 1, padx = c(5, 5), pady = c(5, 5))
    tkgrid(tolerance_ppm_entry, row = 10, column = 2, padx = c(5, 5), pady = c(5, 5))
    tkgrid(commit_preprocessing_button, row = 11, column = 1, columnspan = 3, padx = c(5, 5), pady = c(5, 5))
  }
  
  ##### File type (export)
  file_type_export_choice <- function() {
    # Catch the value from the menu
    file_type_export <- select.list(c("csv","xlsx","xls"), title = "Choose", multiple = FALSE, preselect = "csv")
    # Raise the focus on the main window
    tkraise(window)
    # Default
    if (file_type_export == "") {
      file_type_export <- "csv"
    }
    if (file_type_export == "xls" || file_type_export == "xlsx") {
      # Try to install the XLConnect (it will fail if Java is not installed)
      Java_is_installed <- FALSE
      try({
        install_and_load_required_packages("XLConnect")
        Java_is_installed <- TRUE
      }, silent = TRUE)
      # If it didn't install successfully, set to CSV
      if (Java_is_installed == FALSE) {
        tkmessageBox(title = "Java not installed", message = "Java is not installed, therefore the package XLConnect cannot be installed and loaded.\nThe output format is switched back to CSV", icon = "warning")
        file_type_export <- "csv"
      }
    }
    # Escape the function
    file_type_export <<- file_type_export
    # Set the value of the displaying label
    file_type_export_value_label <- tklabel(window, text = file_type_export, font = label_font, bg = "white", width = 20)
    tkgrid(file_type_export_value_label, row = 8, column = 6, padx = c(10, 10), pady = c(10, 10))
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### File name (export)
  set_file_name <- function() {
    # Retrieve the peaklist file name from the entry...
    filename <- tclvalue(file_name)
    # Add the date and time to the filename
    current_date <- unlist(strsplit(as.character(Sys.time()), " "))[1]
    current_date_split <- unlist(strsplit(current_date, "-"))
    current_time <- unlist(strsplit(as.character(Sys.time()), " "))[2]
    current_time_split <- unlist(strsplit(current_time, ":"))
    final_date <- ""
    for (x in 1:length(current_date_split)) {
      final_date <- paste(final_date, current_date_split[x], sep="")
    }
    final_time <- ""
    for (x in 1:length(current_time_split)) {
      final_time <- paste(final_time, current_time_split[x], sep="")
    }
    final_date_time <- paste(final_date, final_time, sep = "_")
    filename <- paste(filename, " (", final_date_time, ")", sep = "")
    # Create a copy for the subfolder name (for the spectral files)
    filename_subfolder <- filename
    # Add the extension if it is not present in the filename
    if (file_type_export == "csv") {
      if (length(grep(".csv", filename, fixed = TRUE)) == 1) {
        filename <- filename
      }    else {filename <- paste(filename, ".csv", sep="")}
    }
    if (file_type_export == "xlsx") {
      if (length(grep(".xlsx", filename, fixed = TRUE)) == 1) {
        filename <- filename
      }    else {filename <- paste(filename, ".xlsx", sep="")}
    }
    if (file_type_export == "xls") {
      if (length(grep(".xls", filename, fixed = TRUE)) == 1) {
        filename <- filename
      }    else {filename <- paste(filename, ".xls", sep="")}
    }
    # Set the value for displaying purposes
    filename_value <- filename
    #### Exit the function and put the variable into the R workspace
    filename <<- filename
    filename_value <<- filename_value
    filename_subfolder <<- filename_subfolder
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Dump parameters
  dump_parameters <- function() {
    parameter_vector <- c(spectra_format_value, peak_picking_algorithm_value, SNR_value, signals_to_take_value, peak_deisotoping_enveloping_value, low_intensity_peak_removal_threshold_percent_value, low_intensity_peak_removal_threshold_method, peak_filtering_threshold_percentage_value, peak_filtering_mode_value, average_replicates_value, mass_range_value, tof_mode_value, tolerance_ppm_value, transform_data_value, smoothing_value, baseline_subtraction_value, normalization_value, spectral_alignment_value, filepath_import_value, output_folder, filename_value, file_type_export, signals_avg_and_sd_value)
    names(parameter_vector) <- c("Spectra format", "Peak picking algorithm", "S/N", "Most intense signals to take", "Peak deisotoping / enveloping", "Low-intensity peak removal threshold percentage", "Low-intensity peak removal method", "Peak filtering threshold percentage", "Peak filtering mode", "Average replicates", "Mass range", "TOF mode", "Tolerance (in ppm)", "Data transformation", "Smoothing", "Baseline subtraction", "Normalization", "Spectral alignment", "Spectra folder", "Output folder", "File name", "File type", "Signal number statistics")
    parameters_matrix <- cbind(parameter_vector)
    rownames(parameters_matrix) <- names(parameter_vector)
    colnames(parameters_matrix) <- "Parameter value"
    parameters_matrix <<- parameters_matrix
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Samples
  select_samples_function <- function() {
    setwd(getwd())
    ########## Prompt if a folder has to be selected or a single file
    # Catch the value from the popping out menu
    spectra_input_type <- select.list(c("file","folder"), title = "Folder or file?", multiple = FALSE, preselect = "folder")
    if (spectra_input_type == "") {
      spectra_input_type <- "file"
    }
    if (spectra_input_type == "folder") {
      filepath_import_select <- tkmessageBox(title = "Samples", message = "Select the folder for the spectra to be imported.\nIf there are spectral (imzML, txt, csv) files in the folder, each one of them should contain spectra from the same class.\nIf there are subfolders, each one of them should contain spectra from the same class as imzML/txt/csv files", icon = "info")
      filepath_import <- tclvalue(tkchooseDirectory())
      if (!nchar(filepath_import)) {
        tkmessageBox(message = "No folder selected")
      } else {
        tkmessageBox(message = paste("The sample spectra will be read from:", filepath_import))
      }
    } else if (spectra_input_type == "file") {
      filepath_import_select <- tkmessageBox(title = "Samples", message = "Select the file for the spectra to be imported", icon = "info")
      # Filter openable files according to the format
      if (spectra_format == "imzML") {
        filepath_import <- tclvalue(tkgetOpenFile(filetypes = "{{imzML files} {.imzML}}"))
      } else if (spectra_format == "txt") {
        filepath_import <- tclvalue(tkgetOpenFile(filetypes = "{{TXT files} {.txt}}"))
      } else if (spectra_format == "csv") {
        filepath_import <- tclvalue(tkgetOpenFile(filetypes = "{{CSV files} {.csv}}"))
      }
      if (!nchar(filepath_import)) {
        tkmessageBox(message = "No file selected")
      } else {
        tkmessageBox(message = paste("The sample spectra will be read from:", filepath_import))
      }
    }
    # Set the value for displaying purposes
    filepath_import_value <- filepath_import
    # Exit the function and put the variable into the R workspace
    filepath_import <<- filepath_import
    filepath_import_value <<- filepath_import_value
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Output
  browse_output_function <- function() {
    output_folder <- tclvalue(tkchooseDirectory())
    if (!nchar(output_folder)) {
      # Get the output folder from the default working directory
      output_folder <- getwd()
    }
    tkmessageBox(message = paste("Every file will be saved in", output_folder))
    setwd(output_folder)
    # Exit the function and put the variable into the R workspace
    output_folder <<- output_folder
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Exit
  end_session_function <- function () {
    q(save="no")
  }
  
  ##### Peak picking algorithm
  peak_picking_algorithm_choice <- function() {
    # Catch the value from the menu
    peak_picking_algorithm <- select.list(c("SuperSmoother", "MAD"), title = "Choose", multiple = FALSE, preselect = "SuperSmoother")
    # Default
    if (peak_picking_algorithm == "") {
      peak_picking_algorithm <- "SuperSmoother"
    }
    # Set the value of the displaying label
    peak_picking_algorithm_value <- peak_picking_algorithm
    if (peak_picking_algorithm_value == "MAD") {
      peak_picking_algorithm_value <- "Median\nAbsolute Deviation"
    } else if (peak_picking_algorithm_value == "SuperSmoother") {
      peak_picking_algorithm_value <- "Super Smoother"
    }
    peak_picking_algorithm_value_label <- tklabel(window, text = peak_picking_algorithm_value, font = label_font, bg = "white", width = 20, height = 2)
    tkgrid(peak_picking_algorithm_value_label, row = 2, column = 4, padx = c(10, 10), pady = c(10, 10))
    # Escape the function
    peak_picking_algorithm <<- peak_picking_algorithm
    peak_picking_algorithm_value <<- peak_picking_algorithm_value
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Peaks deisotoping or enveloping
  peak_deisotoping_enveloping_choice <- function() {
    # Catch the value from the menu
    peak_deisotoping_enveloping <- select.list(c("Peak Deisotoping","Peak Enveloping", "None"), title = "Peak Deisotoping/Enveloping", multiple = FALSE, preselect = "Peak Deisotoping")
    # Default
    if (peak_deisotoping_enveloping == "") {
      peak_deisotoping_enveloping <- "Peak Deisotoping"
    }
    if (peak_deisotoping_enveloping == "Peak Deisotoping") {
      peak_deisotoping <- TRUE
      peak_enveloping <- FALSE
    } else if (peak_deisotoping_enveloping == "Peak Enveloping") {
      peak_deisotoping <- FALSE
      peak_enveloping <- TRUE
    } else if (peak_deisotoping_enveloping == "None") {
      peak_deisotoping <- FALSE
      peak_enveloping <- FALSE
    }
    # Set the value of the displaying label
    peak_deisotoping_enveloping_value <- peak_deisotoping_enveloping
    peak_deisotoping_enveloping_value_label <- tklabel(window, text = peak_deisotoping_enveloping_value, font = label_font, bg = "white", width = 20)
    tkgrid(peak_deisotoping_enveloping_value_label, row = 3, column = 5, padx = c(10, 10), pady = c(10, 10))
    # Escape the function
    peak_deisotoping <<- peak_deisotoping
    peak_enveloping <<- peak_enveloping
    peak_deisotoping_enveloping_value <<- peak_deisotoping_enveloping_value
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Multicore processing
  allow_parallelization_choice <- function() {
    ##### Messagebox
    tkmessageBox(title = "Parallel processing is resource hungry", message = "Parallel processing is resource hungry.\nBy activating it, the computation becomes faster, but the program will eat a lot of RAM, possibly causing your computer to freeze. If you want to play safe, do not enable it", icon = "warning")
    # Catch the value from the menu
    allow_parallelization <- select.list(c("YES","NO"), title = "Parallelization", multiple = FALSE, preselect = "NO")
    # Default
    if (allow_parallelization == "YES") {
      if (Sys.info()[1] == "Windows") {
        allow_parallelization <- "foreach"
      } else {
        allow_parallelization <- "lapply"
      }
    }
    if (allow_parallelization == "NO" || allow_parallelization == "") {
      allow_parallelization <- FALSE
    }
    # Set the value of the displaying label
    if (allow_parallelization == "foreach" || allow_parallelization == "lapply") {
      allow_parallelization_value <- "YES"
    } else {
      allow_parallelization_value <- "NO"
    }
    allow_parallelization_value_label <- tklabel(window, text = allow_parallelization_value, font = label_font, bg = "white", width = 20)
    tkgrid(allow_parallelization_value_label, row = 7, column = 4)
    # Escape the function
    allow_parallelization <<- allow_parallelization
    allow_parallelization_value <<- allow_parallelization_value
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Average the replicates
  average_replicates_choice <- function() {
    # Catch the value from the menu
    average_replicates <- select.list(c("YES","NO"), title="Choose", multiple = FALSE, preselect = "YES")
    # Default
    if (average_replicates == "YES" || average_replicates == "") {
      average_replicates <- TRUE
      # Select the averaging mode
      averaging_method_input <- select.list(c("Average Spectrum", "Skyline Spectrum", "Median Spectrum"), title = "Choose the averaging method", preselect = "Average Spectrum")
      # Determine the value
      if (averaging_method_input == "Average Spectrum" || averaging_method_input == "") {
        averaging_method <- "mean"
      } else if (averaging_method_input == "Skyline Spectrum") {
        averaging_method <- "sum"
      } else if (averaging_method_input == "Median Spectrum") {
        averaging_method <- "median"
      }
    }
    if (average_replicates == "NO") {
      average_replicates <- FALSE
    }
    # Set the value of the displaying label
    if (average_replicates == TRUE) {
      average_replicates_value <- paste0("YES\n", "(", averaging_method_input, ")") 
    } else {
      average_replicates_value <- "NO"
    }
    average_replicates_value_label <- tklabel(window, text = average_replicates_value, font = label_font, bg = "white", width = 20, height = 2)
    tkgrid(average_replicates_value_label, row = 7, column = 2, padx = c(10, 10), pady = c(10, 10))
    # Escape the function
    average_replicates <<- average_replicates
    averaging_method <<- averaging_method
    average_replicates_value <<- average_replicates_value
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Peak filtering mode
  peak_filtering_mode_choice <- function() {
    # Catch the value from the menu
    peak_filtering_mode <- select.list(c("whole dataset","class-wise"), title = "Peak filtering mode", multiple = FALSE, preselect = "whole dataset")
    # Default
    if (peak_filtering_mode == "") {
      peak_filtering_mode <- "whole dataset"
    }
    # Set the value of the displaying label
    peak_filtering_mode_value <- peak_filtering_mode
    peak_filtering_mode_value_label <- tklabel(window, text = peak_filtering_mode, font = label_font, bg = "white", width = 20)
    tkgrid(peak_filtering_mode_value_label, row = 5, column = 5, padx = c(10, 10), pady = c(10, 10))
    # Escape the function
    peak_filtering_mode <<- peak_filtering_mode
    peak_filtering_mode_value <<- peak_filtering_mode_value
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Low intensity peaks removal Method
  low_intensity_peak_removal_threshold_method_choice <- function() {
    # Catch the value from the menu
    low_intensity_peak_removal_threshold_method <- select.list(c("whole","element-wise"), title="Choose", multiple = FALSE, preselect = "element-wise")
    # Default
    if (low_intensity_peak_removal_threshold_method == "") {
      low_intensity_peak_removal_threshold_method <- "element-wise"
    }
    # Set the value of the displaying label
    low_intensity_peak_removal_threshold_method_value <- low_intensity_peak_removal_threshold_method
    if (low_intensity_peak_removal_threshold_method_value == "whole") {
      low_intensity_peak_removal_threshold_method <- "whole"
    }
    low_intensity_peak_removal_threshold_method_value_label <- tklabel(window, text = low_intensity_peak_removal_threshold_method_value, font = label_font, bg = "white", width = 20)
    tkgrid(low_intensity_peak_removal_threshold_method_value_label, row = 4, column = 5, padx = c(10, 10), pady = c(10, 10))
    # Escape the function
    low_intensity_peak_removal_threshold_method <<- low_intensity_peak_removal_threshold_method
    low_intensity_peak_removal_threshold_method_value <<- low_intensity_peak_removal_threshold_method_value
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### File format
  spectra_format_choice <- function() {
    # Catch the value from the menu
    spectra_format <- select.list(c("imzML", "Xmass", "TXT", "CSV", "MSD"), title = "Spectra format", preselect = "Xmass")
    # Default
    if (spectra_format == "Xmass") {
      spectra_format <- "fid"
      spectra_format_value <- "Xmass"
    } else if (spectra_format == "" || spectra_format == "imzML") {
      spectra_format <- "imzML"
      spectra_format_value <- "imzML"
    } else if (spectra_format == "TXT") {
      spectra_format <- "txt"
      spectra_format_value <- "TXT"
    } else if (spectra_format == "CSV") {
      spectra_format <- "csv"
      spectra_format_value <- "CSV"
    } else if (spectra_format == "MSD") {
      spectra_format <- "msd"
      spectra_format_value <- "MSD"
    }
    # Advisory messages
    if (spectra_format == "txt") {
      tkmessageBox(title = "TXT format", message = "The TXT file should have two columns: the m/z values and the intensity values, separated by a tab and without any header", icon = "info")
    } else if (spectra_format == "csv") {
      tkmessageBox(title = "CSV format", message = "The CSV file should have two columns: the m/z values and the intensity values, separated by a comma and with a header", icon = "info")
    }
    # Escape the function
    spectra_format <<- spectra_format
    spectra_format_value <<- spectra_format_value
    # Set the value of the displaying label
    spectra_format_value_label <- tklabel(window, text = spectra_format_value, font = label_font, bg = "white", width = 20)
    tkgrid(spectra_format_value_label, row = 2, column = 2, padx = c(10, 10), pady = c(10, 10))
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Import the spectra
  import_spectra_function <- function() {
    ##### Run only if the spectra path has been set!
    if (!is.null(filepath_import)) {
      ## Initialization
      spectra <- NULL
      ### Put all the import block under the try() statement, so that if there are blocking errors (such as no files), the spectra variable remains NULL.
      try({
        # Progress bar
        import_progress_bar <- tkProgressBar(title = "Importing and preprocessing spectra...", label = "", min = 0, max = 1, initial = 0, width = 300)
        setTkProgressBar(import_progress_bar, value = 0, title = NULL, label = "0 %")
        ###### Get the values
        # Generate the list of spectra
        if (spectra_format == "fid" || spectra_format == "txt" || spectra_format == "csv" || spectra_format == "msd") {
          ### Load the spectra
          setTkProgressBar(import_progress_bar, value = 0.25, title = NULL, label = "25 %")
          spectra <- import_spectra(filepath = filepath_import, spectra_format = spectra_format, mass_range = mass_range, allow_parallelization = allow_parallelization, spectral_names = "name", replace_sample_name_field = FALSE, remove_empty_spectra = TRUE)
          # Preprocessing
          setTkProgressBar(import_progress_bar, value = 0.50, title = NULL, label = "50 %")
          spectra <- preprocess_spectra(spectra, tof_mode = tof_mode, preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_algorithm_parameter = baseline_subtraction_algorithm_parameter, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, spectral_alignment_algorithm = NULL, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of), allow_parallelization = allow_parallelization, tolerance_ppm = tolerance_ppm)
          setTkProgressBar(import_progress_bar, value = 0.75, title = NULL, label = "75 %")
        } else if (spectra_format == "imzML") {
          # List all the imzML files (if the path is not already an imzML file)
          if (length(grep(".imzML", filepath_import, fixed = TRUE)) <= 0) {
            imzml_files <- read_spectra_files(filepath_import, spectra_format = spectra_format, full_path = TRUE)
          } else {
            imzml_files <- filepath_import
          }
          # Generate the spectra list
          spectra <- NULL
          ### Load the spectra
          if (!is.null(mass_range)) {
            # Read and import one imzML file at a time
            if (length(imzml_files) > 0) {
              setTkProgressBar(import_progress_bar, value = 0.50, title = NULL, label = "50 %")
              for (imzml in 1:length(imzml_files)) {
                # Read and import the imzML file
                spectra_imzml <- importImzMl(imzml_files[imzml], massRange = mass_range)
                # Preprocessing
                spectra_imzml <- preprocess_spectra(spectra_imzml, tof_mode = tof_mode, preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_algorithm_parameter = baseline_subtraction_algorithm_parameter, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, spectral_alignment_algorithm = NULL, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of), allow_parallelization = allow_parallelization, tolerance_ppm = tolerance_ppm)
                # Average the replicates (one AVG spectrum for each imzML file)
                if (average_replicates == TRUE) {
                  spectra_imzml <- averageMassSpectra(spectra_imzml, method = averaging_method)
                  # Preprocessing AVG
                  spectra_imzml <- preprocess_spectra(spectra_imzml, tof_mode = tof_mode, preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_algorithm_parameter = baseline_subtraction_algorithm_parameter, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, spectral_alignment_algorithm = NULL, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of), allow_parallelization = allow_parallelization, tolerance_ppm = tolerance_ppm)
                }
                # Append it to the final list of spectra
                if (is.null(spectra)) {
                  spectra <- spectra_imzml
                } else {
                  spectra <- append(spectra, spectra_imzml)
                }
              }
              setTkProgressBar(import_progress_bar, value = 0.75, title = NULL, label = "75 %")
            }
          } else {
            # Read and import one imzML file at a time
            if (length(imzml_files) > 0) {
              setTkProgressBar(import_progress_bar, value = 0.50, title = NULL, label = "50 %")
              for (imzml in 1:length(imzml_files)) {
                # Read and import the imzML file
                spectra_imzml <- importImzMl(imzml_files[imzml])
                # Preprocessing
                spectra_imzml <- preprocess_spectra(spectra_imzml, tof_mode = tof_mode, preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_algorithm_parameter = baseline_subtraction_algorithm_parameter, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, spectral_alignment_algorithm = NULL, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of), allow_parallelization = allow_parallelization, tolerance_ppm = tolerance_ppm)
                # Average the replicates (one AVG spectrum for each imzML file)
                if (average_replicates == TRUE) {
                  spectra_imzml <- averageMassSpectra(spectra_imzml, method = averaging_method)
                  # Preprocessing AVG
                  spectra_imzml <- preprocess_spectra(spectra_imzml, tof_mode = tof_mode, preprocessing_parameters = list(mass_range = NULL, transformation_algorithm = transform_data_algorithm, smoothing_algorithm = smoothing_algorithm, smoothing_strength = smoothing_strength, baseline_subtraction_algorithm = baseline_subtraction_algorithm, baseline_subtraction_algorithm_parameter = baseline_subtraction_algorithm_parameter, normalization_algorithm = normalization_algorithm, normalization_mass_range = normalization_mass_range, spectral_alignment_algorithm = NULL, preprocess_spectra_in_packages_of = preprocess_spectra_in_packages_of), allow_parallelization = allow_parallelization, tolerance_ppm = tolerance_ppm)
                }
                # Append it to the final list of spectra
                spectra <- append(spectra, spectra_imzml)
              }
              setTkProgressBar(import_progress_bar, value = 0.75, title = NULL, label = "75 %")
            }
          }
          # Replace sample name
          spectra <- replace_sample_name_list(spectra, spectra_format = spectra_format, type = "name", replace_sample_name_field = FALSE)
        }
        ##### Alignment of the imported spectra
        if (!is.null(spectral_alignment_algorithm)) {
          spectral_alignment_performed <- FALSE
          try({
            spectra_names <- names(spectra)
            spectra <- align_spectra(spectra, spectral_alignment_algorithm = spectral_alignment_algorithm, spectral_alignment_reference = spectral_alignment_reference, tof_mode = tof_mode, deisotope_peaklist = FALSE)
            names(spectra) <- spectra_names
            spectral_alignment_performed <- TRUE
          }, silent = TRUE)
          ### Spectral alignment messagebox
          if (spectral_alignment_performed == TRUE) {
            print("The spectral alignment has been succesfully performed!")
          } else {
            print("The spectral alignment could not be performed!")
            #tkmessageBox(title = "Spectral alignment not possible", message = "The spectral alignment could not be performed!", icon = "warning")
          }
        }
        ##### Retrieve the class list
        # List the directories in the filepath_import folder
        folder_list <- list.dirs(filepath_import, full.names = FALSE, recursive = FALSE)
        # If there are only imzML files
        if (length(folder_list) == 0) {
          # Each imzML file is a class
          class_list <- read_spectra_files(filepath_import, spectra_format = spectra_format, full_path = TRUE)
          if (length(class_list) == 0) {
            class_list <- filepath_import
          }
        } else if ((length(folder_list) == 1 && folder_list != "") || (length(folder_list) >= 1)) {
          class_list <- folder_list
        }
        setTkProgressBar(import_progress_bar, value = 1, title = NULL, label = "100 %")
        close(import_progress_bar)
        # Exit the function and put the variable into the R workspace
        spectra <<- spectra
        class_list <<- class_list
      }, silent = TRUE)
      if (is.null(spectra)) {
        try(close(import_progress_bar), silent = TRUE)
        tkmessageBox(title = "No spectral files", message = "There are no spectral files in the selected folder!\n\nTry to select another folder or another format!", icon = "warning")
      } else {
        tkmessageBox(title = "Import successful", message = "The spectra have been successfully imported and preprocessed", icon = "info")
      }
    } else {
      ### Messagebox
      tkmessageBox(title = "Import not possible", message = "No spectra files or folder have been selected!", icon = "warning")
      import_successful <<- FALSE
    }
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Peak picking function
  peak_picking_function <- function() {
    ############ Do not run if the spectra have not been imported
    if (!is.null(spectra)) {
      # Progress bar
      peak_picking_progress_bar <- tkProgressBar(title = "Peak picking...", label = "", min = 0, max = 1, initial = 0, width = 300)
      setTkProgressBar(peak_picking_progress_bar, value = 0, title = NULL, label = "0 %")
      ###### Get the values
      ## Signals to take in most intense peaks
      signals_to_take <- tclvalue(signals_to_take)
      signals_to_take <- as.integer(signals_to_take)
      signals_to_take_value <- as.character(signals_to_take)
      ## SNR
      SNR <- tclvalue(SNR)
      SNR <- as.numeric(SNR)
      SNR_value <- as.character(SNR)
      setTkProgressBar(peak_picking_progress_bar, value = 0.25, title = NULL, label = "25 %")
      peaks <- peak_picking(spectra, peak_picking_algorithm = peak_picking_algorithm, SNR = SNR, tof_mode = tof_mode, allow_parallelization = allow_parallelization, deisotope_peaklist = peak_deisotoping, envelope_peaklist = peak_enveloping, signals_to_take = signals_to_take)
      ## Peaks filtering threshold
      peak_filtering_threshold_percentage <- tclvalue(peak_filtering_threshold_percentage)
      peak_filtering_threshold_percentage <- as.numeric(peak_filtering_threshold_percentage)
      peak_filtering_threshold_percentage_value <- as.character(peak_filtering_threshold_percentage)
      ## Low intensity threshold
      low_intensity_peak_removal_threshold_percent <- tclvalue(low_intensity_peak_removal_threshold_percent)
      low_intensity_peak_removal_threshold_percent <- as.numeric(low_intensity_peak_removal_threshold_percent)
      low_intensity_peak_removal_threshold_percent_value <- as.character(low_intensity_peak_removal_threshold_percent)
      ##### Peak alignment
      setTkProgressBar(peak_picking_progress_bar, value = 0.50, title = NULL, label = "50 %")
      if (peak_filtering_mode == "class-wise") {
        if (isMassPeaksList(peaks)) {
          peaks_class <- replace_class_name(peaks, class_list = class_list, class_in_file_path = TRUE, class_in_file_name = FALSE, spectra_format = spectra_format)
          class_vector_for_peak_filtering <- vector()
          for (p in 1:length(peaks_class)) {
            class_vector_for_peak_filtering <- append(class_vector_for_peak_filtering, peaks_class[[p]]@metaData$file)
          }
          peaks <- align_and_filter_peaks(peaks, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peak_filtering_threshold_percentage, class_vector_for_peak_filtering = class_vector_for_peak_filtering, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, reference_peaklist = NULL, spectra = spectra, alignment_iterations = 5, allow_parallelization = allow_parallelization)
        }
      } else if (peak_filtering_mode == "whole dataset") {
        peaks <- align_and_filter_peaks(peaks, peak_picking_algorithm = peak_picking_algorithm, tof_mode = tof_mode, peak_filtering_frequency_threshold_percent = peak_filtering_threshold_percentage, class_vector_for_peak_filtering = NULL, low_intensity_peak_removal_threshold_percent = low_intensity_peak_removal_threshold_percent, low_intensity_peak_removal_threshold_method = low_intensity_peak_removal_threshold_method, reference_peaklist = NULL, spectra = spectra, alignment_iterations = 5, allow_parallelization = allow_parallelization)
      }
      setTkProgressBar(peak_picking_progress_bar, value = 1, title = NULL, label = "100 %")
      close(peak_picking_progress_bar)
      # Exit the function and put the variable into the R workspace
      peaks <<- peaks
      signals_to_take_value <<- signals_to_take_value
      SNR_value <<- SNR_value
      peak_filtering_threshold_percentage_value <<- peak_filtering_threshold_percentage_value
      low_intensity_peak_removal_threshold_percent_value <<- low_intensity_peak_removal_threshold_percent_value
      ### Messagebox
      tkmessageBox(title = "Peak picking successful", message = "The peak picking process has been successfully performed", icon = "info")
    } else if (is.null(spectra)) {
      ### Messagebox
      tkmessageBox(title = "Spectra not imported", message = "The spectra have not been imported yet.\nImport them before performing the peak picking", icon = "warning")
    }
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Output the average number of signals with the SD
  signals_avg_and_sd_function <- function() {
    ############ Do not run if the peaks have not been picked
    if (!is.null(peaks)) {
      # Generate the vector recording the number of signals
      number_of_signals_vector <- numeric()
      if (isMassPeaksList(peaks)) {
        for (p in 1:length(peaks)) {
          number_of_signals_vector <- append(number_of_signals_vector, length(peaks[[p]]@mass))
        }
        # Compute the mean and the standard deviation
        mean_signal_number <- mean(number_of_signals_vector, na.rm = TRUE)
        median_signal_number <- median(number_of_signals_vector, na.rm = TRUE)
        sd_signal_number <- sd(number_of_signals_vector, na.rm = TRUE)
        cv_signal_number <- (sd_signal_number / mean_signal_number) * 100
        ### Messagebox
        # Message
        message_avg_sd <- paste("The mean number of signals in the spectral dataset is:", mean_signal_number, ",\n\nthe standard deviation is:", sd_signal_number, ",\n\nthe coefficient of variation is:", cv_signal_number, "%")
        tkmessageBox(title = "Mean and SD of the number of signals", message = message_avg_sd, icon = "info")
        signals_avg_and_sd_value <- paste0("Mean: ", mean_signal_number, "\nMedian: ", median_signal_number, "\nStandard Deviation: ", sd_signal_number, "\nCoefficient of Variation: ", cv_signal_number)
      } else if (isMassPeaks(peaks)) {
        number_of_signals <- length(peaks@mass)
        # Message
        message_avg_sd <- paste("The number of signals in the spectrum is:", number_of_signals)
        tkmessageBox(title = "Number of signals", message = message_avg_sd, icon = "info")
        signals_avg_and_sd_value <- paste0("Number of signals: ", number_of_signals)
      }
      signals_avg_and_sd_value <<- signals_avg_and_sd_value
      # Generate a label to see the output
      signals_avg_and_sd_value_label <- tklabel(window, text = signals_avg_and_sd_value, font = label_font, bg = "white", width = 40, height = 5)
      tkgrid(signals_avg_and_sd_value_label, row = 10, column = 4, padx = c(10, 10), pady = c(10, 10),columnspan = 2)
    } else if (is.null(peaks)) {
      ### Messagebox
      tkmessageBox(title = "Something is wrong", message = "Some elements are needed to perform this operation: make sure that the peak picking process has been performed", icon = "warning")
    }
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Run the Peaklist Export function
  run_peaklist_export_function <- function() {
    setwd(output_folder)
    # Get the filename from the entry
    set_file_name()
    ######## Run only if all the elements needed are there
    if (!is.null(spectra) && !is.null(peaks)) {
      # Progress bar
      program_progress_bar <- tkProgressBar(title = "Computing...", label = "", min = 0, max = 1, initial = 0, width = 300)
      setTkProgressBar(program_progress_bar, value = 0, title = NULL, label = "0 %")
      setTkProgressBar(program_progress_bar, value = 0.50, title = NULL, label = "50 %")
      # Generate the signal matrix
      peaklist <- intensityMatrix(peaks, spectra)
      setTkProgressBar(program_progress_bar, value = 0.75, title = NULL, label = "75 %")
      # Add the sample and the class to the matrix
      peaklist <- matrix_add_class_and_sample(peaklist, peaks = peaks, class_list = class_list, spectra_format = spectra_format, sample_output = TRUE, class_output = TRUE, row_labels = NULL)
      setTkProgressBar(program_progress_bar, value = 0.90, title = NULL, label = "90 %")
      # Save the files (CSV)
      if (file_type_export == "csv") {
        
        write.csv(peaklist, file = filename, row.names = FALSE)
        dump_parameters()
        write.csv(parameters_matrix, file = paste0(filename_subfolder, " - Parameters.", file_type_export), row.names = TRUE, col.names = TRUE)
      } else if (file_type_export == "xlsx" || file_type_export == "xls") {
        # Save the files (Excel)
        peaklist <- as.data.frame(peaklist)
        # Generate unique row names
        unique_row_names <- make.names(rownames(peaklist), unique = TRUE)
        rownames(peaklist) <- unique_row_names
        # Export
        writeWorksheetToFile(file = filename, data = peaklist, sheet = "Peaklist", clearSheets = TRUE, header = TRUE)
        #write.xlsx(x = peaklist, file = filename, sheetName="Peaklist", row.names = FALSE)
        dump_parameters()
        writeWorksheetToFile(file = paste0(filename_subfolder, " - Parameters.", file_type_export), data = parameters_matrix, sheet = "Parameters", clearSheets = TRUE, header = TRUE)
      }
      setTkProgressBar(program_progress_bar, value = 1, title = NULL, label = "100 %")
      close(program_progress_bar)
      ### Messagebox
      tkmessageBox(title = "Done!", message = "The peaklist file has been dumped!", icon = "info")
    } else if (is.null(spectra) || is.null(peaks)) {
      ### Messagebox
      tkmessageBox(title = "Something is wrong", message = "Some elements are needed to perform this operation: make sure that the spectra have been imported and the peak picking process has been performed", icon = "warning")
    }
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Dump the spectral files
  dump_spectra_files_function <- function() {
    ### Run only if there are spectra
    if (!is.null(spectra) && !is.null(peaks)) {
      # Get the filename from the entry (filename_subfolder)
      # Go to the working directory and create a folder named 'Spectra files'
      spectra_files_subfolder <- file.path(output_folder, paste(filename_subfolder, "- Spectral files"))
      dir.create(spectra_files_subfolder)
      setwd(spectra_files_subfolder)
      # Replace the sample path with the sample name in the metadata
      spectra <- replace_sample_name(spectra, spectra_format = spectra_format, allow_parallelization = allow_parallelization)
      # Get the names of the spectra and generate a vector of names
      spectra_name_vector <- character()
      if (isMassSpectrumList(spectra)) {
        if (is.null(names(spectra))) {
          for (s in 1:length(spectra)) {
            spectra_name_vector <- append(spectra_name_vector, spectra[[s]]@metaData$file[1])
          }
        } else {
          spectra_name_vector <- names(spectra)
        }
      } else if (isMassSpectrum(spectra)) {
        if (is.null(names(spectra))) {
          spectra_name_vector <- spectra@metaData$file[1]
        } else {
          spectra_name_vector <- names(spectra)
        }
      }
      # If there are already unique names, leave the spectra_name_vector as it is...
      if (length(spectra_name_vector) == length(unique(spectra_name_vector))) {
        spectra_name_vector <- spectra_name_vector
      } else {
        # Otherwise, generate unique names...
        spectra_name_vector <- make.names(spectra_name_vector, unique = TRUE)
      }
      ### Dump the spectal files
      # Choose the file format
      spectra_output_format <- select.list(c("MSD", "TXT"), preselect = "MSD", multiple = FALSE, title = "Select the spectra format")
      if (spectra_output_format == "") {
        spectra_output_format <- "MSD"
      }
      # MSD
      if (spectra_output_format == "MSD") {
        if (isMassSpectrumList(spectra)) {
          if (isMassPeaksList(peaks) && length(peaks) == length(spectra)) {
            for (s in 1:length(spectra)) {
              exportMsd(spectra[[s]], file = paste(spectra_name_vector[s], ".msd", sep = ""), force = TRUE, peaks = peaks[[s]])
            }
          } else {
            for (s in 1:length(spectra)) {
              exportMsd(spectra[[s]], file = paste(spectra_name_vector[s], ".msd", sep = ""), force = TRUE)
            }
          }
        } else if (isMassSpectrum(spectra)) {
          if (isMassPeaks(peaks)) {
            exportMsd(spectra, file = paste(spectra_name_vector, ".msd", sep = ""), force = TRUE, peaks = peaks)
          } else {
            exportMsd(spectra, file = paste(spectra_name_vector, ".msd", sep = ""), force = TRUE)
          }
        }
      } else if (spectra_output_format == "TXT") {
        if (isMassSpectrumList(spectra)) {
          if (isMassPeaksList(peaks) && length(peaks) == length(spectra)) {
            for (s in 1:length(spectra)) {
              spectra_txt <- matrix(0, ncol = 2, nrow = length(spectra[[s]]@mass))
              peaks_txt <- matrix(0, ncol = 2, nrow = length(peaks[[s]]@mass))
              spectra_txt[, 1] <- cbind(spectra[[s]]@mass)
              spectra_txt[, 2] <- cbind(spectra[[s]]@intensity)
              peaks_txt[, 1] <- cbind(peaks[[s]]@mass)
              peaks_txt[, 2] <- cbind(peaks[[s]]@intensity)
              write.table(spectra_txt, file = paste(spectra_name_vector[s], ".txt", sep = ""), row.names = FALSE, col.names = FALSE)
              write.table(peaks_txt, file = paste(spectra_name_vector[s], " - Peaks.txt", sep = ""), row.names = FALSE, col.names = FALSE)
            }
          } else {
            for (s in 1:length(spectra)) {
              spectra_txt <- matrix(0, ncol = 2, nrow = length(spectra[[s]]@mass))
              spectra_txt[, 1] <- cbind(spectra[[s]]@mass)
              spectra_txt[, 2] <- cbind(spectra[[s]]@intensity)
              write.table(spectra_txt, file = paste(spectra_name_vector[s], ".txt", sep = ""), row.names = FALSE, col.names = FALSE)
            }
          }
        } else if (isMassSpectrum(spectra)) {
          spectra_txt <- matrix(0, ncol = 2, nrow = length(spectra@mass))
          peaks_txt <- matrix(0, ncol = 2, nrow = length(peaks@mass))
          spectra_txt[, 1] <- cbind(spectra@mass)
          spectra_txt[, 2] <- cbind(spectra@intensity)
          peaks_txt[, 1] <- cbind(peaks@mass)
          peaks_txt[, 2] <- cbind(peaks@intensity)
          write.table(spectra_txt, file = paste(spectra_name_vector, ".txt", sep = ""), row.names = FALSE, col.names = FALSE)
          write.table(peaks_txt, file = paste(spectra_name_vector, " - Peaks.txt", sep = ""), row.names = FALSE, col.names = FALSE)
        }
      }
      # Go back to the output folder
      setwd(output_folder)
      ### Messagebox
      tkmessageBox(title = "Spectra files dumped", message = "The spectra files have been succesfully dumped!", icon = "info")
    } else {
      ### Messagebox
      tkmessageBox(title = "Spectra not loaded or Peaks not picked!", message = "No spectra have been imported yet or no peak picking has been performed!", icon = "warning")
    }
    # Raise the focus on the main window
    tkraise(window)
  }
  
  ##### Show info function
  show_info_function <- function() {
    if (Sys.info()[1] == "Linux") {
      system(command = paste("xdg-open", github_wiki_url), intern = FALSE)
    } else if (Sys.info()[1] == "Darwin") {
      system(command = paste("open", github_wiki_url), intern = FALSE)
    } else if (Sys.info()[1] == "Windows") {
      system(command = paste("cmd /c start", github_wiki_url), intern = FALSE)
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  ##################################################################### WINDOW GUI
  
  ########## List of variables, whose values are taken from the entries in the GUI
  SNR <- tclVar("")
  peak_filtering_threshold_percentage <- tclVar("")
  low_intensity_peak_removal_threshold_percent <- tclVar("")
  signals_to_take <- tclVar("")
  file_name <- tclVar("")
  
  
  
  ######################## GUI
  
  ### Get system info (Platform - Release - Version (- Linux Distro))
  system_os = Sys.info()[1]
  os_release = Sys.info()[2]
  os_version = Sys.info()[3]
  
  ### Get the screen resolution
  try({
    # Windows
    if (system_os == "Windows") {
      # Get system info
      screen_info <- system("wmic path Win32_VideoController get VideoModeDescription", intern = TRUE)[2]
      # Get the resolution
      screen_resolution <- unlist(strsplit(screen_info, "x"))
      # Retrieve the values
      screen_height <- as.numeric(screen_resolution[2])
      screen_width <- as.numeric(screen_resolution[1])
    } else if (system_os == "Linux") {
      # Get system info
      screen_info <- system("xdpyinfo -display :0", intern = TRUE)
      # Get the resolution
      screen_resolution <- screen_info[which(screen_info == "screen #0:") + 1]
      screen_resolution <- unlist(strsplit(screen_resolution, "dimensions: ")[1])
      screen_resolution <- unlist(strsplit(screen_resolution, "pixels"))[2]
      # Retrieve the wto dimensions...
      screen_width <- as.numeric(unlist(strsplit(screen_resolution, "x"))[1])
      screen_height <- as.numeric(unlist(strsplit(screen_resolution, "x"))[2])
    }
  }, silent = TRUE)
  
  
  
  ### FONTS
  # Default sizes (determined on a 1680x1050 screen) (in order to make them adjust to the size screen, the screen resolution should be retrieved)
  title_font_size_default <- 18
  other_font_size_default <- 9
  title_font_size <- title_font_size_default
  other_font_size <- other_font_size_default
  
  # Adjust fonts size according to the pixel number
  try({
    # Windows
    if (system_os == "Windows") {
      # Determine the font size according to the resolution
      total_number_of_pixels <- screen_width * screen_height
      # Determine the scaling factor (according to a complex formula)
      scaling_factor_title_font <- as.numeric((0.03611 * total_number_of_pixels) + 9803.1254)
      scaling_factor_other_font <- as.numeric((0.07757 * total_number_of_pixels) + 23529.8386)
      title_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_title_font))
      other_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_other_font))
    } else if (system_os == "Linux") {
      # Linux
      # Determine the font size according to the resolution
      total_number_of_pixels <- screen_width * screen_height
      # Determine the scaling factor (according to a complex formula)
      scaling_factor_title_font <- as.numeric((0.03611 * total_number_of_pixels) + 9803.1254)
      scaling_factor_other_font <- as.numeric((0.07757 * total_number_of_pixels) + 23529.8386)
      title_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_title_font))
      other_font_size <- as.integer(round(total_number_of_pixels / scaling_factor_other_font))
    } else if (system_os == "Darwin") {
      # macOS
      print("Using default font sizes...")
    }
    # Go back to defaults if there are NAs
    if (is.na(title_font_size)) {
      title_font_size <- title_font_size_default
    }
    if (is.na(other_font_size)) {
      other_font_size <- other_font_size_default
    }
  }, silent = TRUE)
  
  # Define the fonts
  # Windows
  if (system_os == "Windows") {
    garamond_title_bold = tkfont.create(family = "Garamond", size = title_font_size, weight = "bold")
    garamond_other_normal = tkfont.create(family = "Garamond", size = other_font_size, weight = "normal")
    arial_title_bold = tkfont.create(family = "Arial", size = title_font_size, weight = "bold")
    arial_other_normal = tkfont.create(family = "Arial", size = other_font_size, weight = "normal")
    trebuchet_title_bold = tkfont.create(family = "Trebuchet MS", size = title_font_size, weight = "bold")
    trebuchet_other_normal = tkfont.create(family = "Trebuchet MS", size = other_font_size, weight = "normal")
    trebuchet_other_bold = tkfont.create(family = "Trebuchet MS", size = other_font_size, weight = "bold")
	calibri_title_bold = tkfont.create(family = "Calibri", size = title_font_size, weight = "bold")
    calibri_other_normal = tkfont.create(family = "Calibri", size = other_font_size, weight = "normal")
    calibri_other_bold = tkfont.create(family = "Calibri", size = other_font_size, weight = "bold")
    # Use them in the GUI
    title_font = calibri_title_bold
    label_font = calibri_other_normal
    entry_font = calibri_other_normal
    button_font = calibri_other_bold
  } else if (system_os == "Linux") {
    #Linux
    # Ubuntu
    if (length(grep("Ubuntu", os_version, ignore.case = TRUE)) > 0) {
      # Define the fonts
      ubuntu_title_bold = tkfont.create(family = "Ubuntu", size = (title_font_size + 2), weight = "bold")
      ubuntu_other_normal = tkfont.create(family = "Ubuntu", size = (other_font_size), weight = "normal")
      ubuntu_other_bold = tkfont.create(family = "Ubuntu", size = (other_font_size), weight = "bold")
      liberation_title_bold = tkfont.create(family = "Liberation Sans", size = title_font_size, weight = "bold")
      liberation_other_normal = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "normal")
      liberation_other_bold = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "bold")
      bitstream_charter_title_bold = tkfont.create(family = "Bitstream Charter", size = title_font_size, weight = "bold")
      bitstream_charter_other_normal = tkfont.create(family = "Bitstream Charter", size = other_font_size, weight = "normal")
      bitstream_charter_other_bold = tkfont.create(family = "Bitstream Charter", size = other_font_size, weight = "bold")
      # Use them in the GUI
      title_font = ubuntu_title_bold
      label_font = ubuntu_other_normal
      entry_font = ubuntu_other_normal
      button_font = ubuntu_other_bold
    } else if (length(grep("Fedora", os_version, ignore.case = TRUE)) > 0) {
      # Fedora
      cantarell_title_bold = tkfont.create(family = "Cantarell", size = title_font_size, weight = "bold")
      cantarell_other_normal = tkfont.create(family = "Cantarell", size = other_font_size, weight = "normal")
      cantarell_other_bold = tkfont.create(family = "Cantarell", size = other_font_size, weight = "bold")
      liberation_title_bold = tkfont.create(family = "Liberation Sans", size = title_font_size, weight = "bold")
      liberation_other_normal = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "normal")
      liberation_other_bold = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "bold")
      # Use them in the GUI
      title_font = cantarell_title_bold
      label_font = cantarell_other_normal
      entry_font = cantarell_other_normal
      button_font = cantarell_other_bold
    } else {
      # Other linux distros
      liberation_title_bold = tkfont.create(family = "Liberation Sans", size = title_font_size, weight = "bold")
      liberation_other_normal = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "normal")
      liberation_other_bold = tkfont.create(family = "Liberation Sans", size = other_font_size, weight = "bold")
      # Use them in the GUI
      title_font = liberation_title_bold
      label_font = liberation_other_normal
      entry_font = liberation_other_normal
      button_font = liberation_other_bold
    }
  } else if (system_os == "Darwin") {
    # macOS
    helvetica_title_bold = tkfont.create(family = "Helvetica", size = title_font_size, weight = "bold")
    helvetica_other_normal = tkfont.create(family = "Helvetica", size = other_font_size, weight = "normal")
    helvetica_other_bold = tkfont.create(family = "Helvetica", size = other_font_size, weight = "bold")
    # Use them in the GUI
    title_font = helvetica_title_bold
    label_font = helvetica_other_normal
    entry_font = helvetica_other_normal
    button_font = helvetica_other_bold
  }
  
  
  
  # The "area" where we will put our input lines (not resizable)
  window <- tktoplevel(bg = "white")
  #tkpack.propagate(window, FALSE)
  tkwm.resizable(window, FALSE, FALSE)
  # Raise the focus on the main window
  tkraise(window)
  tktitle(window) <- "MS PEAKLIST EXPORT"
  #### Browse
  # Title label
  title_label <- tkbutton(window, text = "MS PEAKLIST EXPORT", command = show_info_function, font = title_font, bg = "white", relief = "flat")
  # Library
  select_samples_button <- tkbutton(window, text="BROWSE\nSPECTRA...", command = select_samples_function, font = button_font, bg = "white", width = 20)
  # Output
  browse_output_button <- tkbutton(window, text="BROWSE\nOUTPUT FOLDER...", command = browse_output_function, font = button_font, bg = "white", width = 20)
  #### Entries
  # Peak picking method
  peak_picking_algorithm_entry <- tkbutton(window, text="PEAK PICKING\nALGORITHM", command = peak_picking_algorithm_choice, font = button_font, bg = "white", width = 20)
  # Signals to take
  signals_to_take_label <- tklabel(window, text="Most intense signals\nto take\n(0 = retain all peaks)", font = button_font, bg = "white", width = 20)
  signals_to_take_entry <- tkentry(window, textvariable = signals_to_take, font = entry_font, bg = "white", width = 5, justify = "center")
  tkinsert(signals_to_take_entry, "end", "0")
  # SNR
  SNR_label <- tklabel(window, text="Signal-to-noise\nratio", font = button_font, bg = "white", width = 20)
  SNR_entry <- tkentry(window, textvariable = SNR, font = entry_font, bg = "white", width = 5, justify = "center")
  tkinsert(SNR_entry, "end", "3")
  # Peaks filtering threshold
  peak_filtering_threshold_percentage_label <- tklabel(window, text="Peak filtering\nthreshold\nfrequency percentage", font = button_font, bg = "white", width = 20)
  peak_filtering_threshold_percentage_entry <- tkentry(window, textvariable = peak_filtering_threshold_percentage, font = entry_font, bg = "white", width = 5, justify = "center")
  tkinsert(peak_filtering_threshold_percentage_entry, "end", "5")
  # Peaks filtering mode
  peak_filtering_mode_entry <- tkbutton(window, text="PEAK FILTERING\nMODE", command = peak_filtering_mode_choice, font = button_font, bg = "white", width = 20)
  # Peaks deisotoping
  peak_deisotoping_entry <- tkbutton(window, text="PEAK\nDEISOTOPING\nor\nENVELOPING", command = peak_deisotoping_enveloping_choice, font = button_font, bg = "white", width = 20)
  # Intensity percentage threshold
  low_intensity_peak_removal_threshold_percent_label <- tklabel(window, text="Low-intensity peak\nremoval\npercentage threshold", font = button_font, bg = "white", width = 20)
  low_intensity_peak_removal_threshold_percent_entry <- tkentry(window, textvariable = low_intensity_peak_removal_threshold_percent, font = entry_font, bg = "white", width = 5, justify = "center")
  tkinsert(low_intensity_peak_removal_threshold_percent_entry, "end", "0")
  # Intensity percentage theshold method
  low_intensity_peak_removal_threshold_method_entry <- tkbutton(window, text="INTENSITY\nTHRESHOLD\nMETHOD", command = low_intensity_peak_removal_threshold_method_choice, font = button_font, bg = "white", width = 20)
  # File format
  spectra_format_entry <- tkbutton(window, text="SPECTRA\nFORMAT", command = spectra_format_choice, font = button_font, bg = "white", width = 20)
  # File type export
  file_type_export_entry <- tkbutton(window, text="FILE TYPE\nEXPORT", command = file_type_export_choice, font = button_font, bg = "white", width = 20)
  # Average the replicates
  average_replicates_button <- tkbutton(window, text="AVERAGE\nREPLICATES", command = average_replicates_choice, font = button_font, bg = "white", width = 20)
  # End session
  #end_session_label <- tklabel(window, text="Quit", font = label_font)
  end_session_button <- tkbutton(window, text="QUIT", command = end_session_function, font = button_font, bg = "white", width = 20)
  # Import the spectra
  import_spectra_button <- tkbutton(window, text="IMPORT AND\nPREPROCESS\nSPECTRA...", command = import_spectra_function, font = button_font, bg = "white", width = 20)
  # Import the spectra
  dump_spectra_files_button <- tkbutton(window, text="DUMP SPECTRA\nFILES...", command = dump_spectra_files_function, font = button_font, bg = "white", width = 20)
  # Peak picking
  peak_picking_button <- tkbutton(window, text="PEAK\nPICKING...", command = peak_picking_function, font = button_font, bg = "white", width = 20)
  # Run the Peaklist Export!!
  peaklist_export_button <- tkbutton(window, text="EXPORT\nPEAKLIST...", command = run_peaklist_export_function, font = button_font, bg = "white", width = 20)
  # Average number of signals and standard deviation
  signals_avg_and_sd_button <- tkbutton(window, text="MEAN +/- SD of\nnumber of signals...", command = signals_avg_and_sd_function, font = button_font, bg = "white", width = 20)
  # Multicore
  allow_parallelization_button <- tkbutton(window, text="ALLOW\nPARALLEL\nCOMPUTING", command = allow_parallelization_choice, font = button_font, bg = "white", width = 20)
  # Spectra preprocessing button
  spectra_preprocessing_button <- tkbutton(window, text="SPECTRA\nPREPROCESSING\nPARAMETERS...", command = preprocessing_window_function, font = button_font, bg = "white", width = 20)
  # Set the file name
  set_file_name_label <- tklabel(window, text="<--- Set the file name", font = label_font, bg = "white", width = 20)
  set_file_name_entry <- tkentry(window, textvariable = file_name, font = entry_font, bg = "white", width = 40, justify = "center")
  tkinsert(set_file_name_entry, "end", "Peaklist")
  # Updates
  download_updates_button <- tkbutton(window, text="DOWNLOAD\nUPDATE...", command = download_updates_function, font = button_font, bg = "white", width = 20)
  
  #### Displaying labels
  file_type_export_value_label <- tklabel(window, text = file_type_export, font = label_font, bg = "white", width = 20)
  peak_picking_algorithm_value_label <- tklabel(window, text = peak_picking_algorithm_value, font = label_font, bg = "white", width = 20, height = 2)
  peak_filtering_mode_value_label <- tklabel(window, text = peak_filtering_mode_value, font = label_font, bg = "white", width = 20)
  peak_deisotoping_enveloping_value_label <- tklabel(window, text = peak_deisotoping_enveloping_value, font = label_font, bg = "white", width = 20)
  low_intensity_peak_removal_threshold_method_value_label <- tklabel(window, text = low_intensity_peak_removal_threshold_method_value, font = label_font, bg = "white", width = 20)
  spectra_format_value_label <- tklabel(window, text = spectra_format_value, font = label_font, bg = "white", width = 20)
  allow_parallelization_value_label <- tklabel(window, text = allow_parallelization_value, font = label_font, bg = "white", width = 20)
  average_replicates_value_label <- tklabel(window, text = average_replicates_value, font = label_font, bg = "white", width = 20, height = 2)
  check_for_updates_value_label <- tkbutton(window, text = check_for_updates_value, command = force_check_for_updates_function, font = label_font, bg = "white", width = 20, relief = "flat")
  
  #### Geometry manager
  # Scrollbar
  #window_scrollbar <- tkscrollbar(window, command = function(...)tkyview(window,...))
  # tkgrid
  tkgrid(title_label, row = 1, column = 1, padx = c(20, 20), pady = c(20, 20), columnspan = 4)
  tkgrid(select_samples_button, row = 8, column = 1, padx = c(10, 10), pady = c(10, 10))
  tkgrid(browse_output_button, row = 8, column = 2, padx = c(10, 10), pady = c(10, 10))
  tkgrid(set_file_name_entry, row = 8, column = 3, padx = c(10, 10), pady = c(10, 10))
  tkgrid(set_file_name_label, row = 8, column = 4, padx = c(10, 10), pady = c(10, 10))
  tkgrid(signals_to_take_label, row = 3, column = 2, padx = c(10, 10), pady = c(10, 10))
  tkgrid(signals_to_take_entry, row = 3, column = 3, padx = c(10, 10), pady = c(10, 10))
  tkgrid(SNR_label, row = 2, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(SNR_entry, row = 2, column = 6, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_filtering_threshold_percentage_label, row = 5, column = 2, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_filtering_threshold_percentage_entry, row = 5, column = 3, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_filtering_mode_entry, row = 5, column = 4, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_filtering_mode_value_label, row = 5, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_deisotoping_entry, row = 3, column = 4, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_deisotoping_enveloping_value_label, row = 3, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(low_intensity_peak_removal_threshold_percent_label, row = 4, column = 2, padx = c(10, 10), pady = c(10, 10))
  tkgrid(low_intensity_peak_removal_threshold_percent_entry, row = 4, column = 3, padx = c(10, 10), pady = c(10, 10))
  tkgrid(low_intensity_peak_removal_threshold_method_entry, row = 4, column = 4, padx = c(10, 10), pady = c(10, 10))
  tkgrid(low_intensity_peak_removal_threshold_method_value_label, row = 4, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_picking_algorithm_entry, row = 2, column = 3, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_picking_algorithm_value_label, row = 2, column = 4, padx = c(10, 10), pady = c(10, 10))
  tkgrid(spectra_format_entry, row = 2, column = 1, padx = c(10, 10), pady = c(10, 10))
  tkgrid(spectra_format_value_label, row = 2, column = 2, padx = c(10, 10), pady = c(10, 10))
  tkgrid(file_type_export_entry, row = 8, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(file_type_export_value_label, row = 8, column = 6, padx = c(10, 10), pady = c(10, 10))
  tkgrid(average_replicates_button, row = 7, column = 1, padx = c(10, 10), pady = c(10, 10))
  tkgrid(average_replicates_value_label, row = 7, column = 2, padx = c(10, 10), pady = c(10, 10))
  tkgrid(allow_parallelization_button, row = 7, column = 3, padx = c(10, 10), pady = c(10, 10))
  tkgrid(allow_parallelization_value_label, row = 7, column = 4, padx = c(10, 10), pady = c(10, 10))
  tkgrid(spectra_preprocessing_button, row = 7, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(import_spectra_button, row = 9, column = 1, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peak_picking_button, row = 9, column = 2, padx = c(10, 10), pady = c(10, 10))
  tkgrid(peaklist_export_button, row = 9, column = 3, padx = c(10, 10), pady = c(10, 10))
  tkgrid(signals_avg_and_sd_button, row = 9, column = 4, padx = c(10, 10), pady = c(10, 10))
  tkgrid(dump_spectra_files_button, row = 9, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(end_session_button, row = 9, column = 6, padx = c(10, 10), pady = c(10, 10))
  tkgrid(download_updates_button, row = 1, column = 5, padx = c(10, 10), pady = c(10, 10))
  tkgrid(check_for_updates_value_label, row = 1, column = 6, padx = c(10, 10), pady = c(10, 10))
  
  set_file_name()
  
  
  
  ################################################################################
}





### Call the functions
functions_mass_spectrometry()

### Run the function
ms_peaklist_export()
