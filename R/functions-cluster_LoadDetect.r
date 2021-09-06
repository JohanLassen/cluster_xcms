


# pd must be a tsv file containing following columns:
# | path
load_n_detect_filewise <- function(pd, file_name, data_folder) {


    if (typeof(a) == "character") {
        phenodata <- readr::read_tsv(pd)
    } else {
        phenodata <- pd; rm(pd)
    }

    # Filename is the location of the single file
    phenodata <- phenodata %>%
        dplyr::filter(path == file_name) %>%
        dplyr::select(sample_name, sample_group)

    # Make XCMSnExp object
    loaded_file <- readMSData(files = phenodata$filename[1],
                            pdata = new("NAnnotatedDataFrame", phenodata),
                            mode = "onDisk",
                            msLevel. = 1L)

    name <- gsub(data_folder, "", file_name)
    name <- gsub("/", "_", file_name)

    output_file <- paste0("./peaks_identified/", gsub("[.].*", ".Rdata", name))

    cwp        <- CentWaveParam(peakwidth = c(4, 30),
                                snthresh = 6,
                                ppm = 12,
                                prefilter = c(3, 200),
                                mzdiff = 0.01)

    peaks_file <- findChromPeaks(loaded_file, cwp)

    save(peaks_file, file = output_file)

    rm(peaks_file, cwp, output_file, name, loaded_file, phenodata, )
}

# The data folder is the top branch folder of all the files.
# It is made like this to ensure that any information from subfolder naming are kept in the filenames during further analysis

load_n_detect_all <- function(pd, data_folder) {

    pd <- readr::read_tsv(pd)
    files <- pd$path

    for (i in seq_along(files)) {
        
        file <- files[i]

        # Check if sample already has been loaded and peak detected before
        name <- gsub(data_folder, "", file)
        name <- gsub("/", "_", file)
        output_file <- paste0("./peaks_identified/", gsub("[.].*", ".Rdata", name))
        if (file.exists(output_file)) {
            next
        }

        # If file has not been loaded, do it now.
        cat("Loading and detecting peaks of file", i, "of", length(files), "")
        load_n_detect_filewise(phenodata = pd, file_name = file, data_folder=data_folder)
        }

}
