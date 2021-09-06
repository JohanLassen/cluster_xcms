

### Concatenate files from peak load- and detection-step

concatenate <- function(pd) {
    files <- read_tsv(pd)
    files <- gsub("./data/", "", files$path)
    files <- gsub("/", "_", files)
    files <- gsub("[.].*", ".Rdata", files)


    concatenate_string <- "c("
    for (i in 1:length(files)) {
        cat(i, " of ", length(files), "\n") #progress bar
        name <- paste0("peak", i)
        file <- files[i]

        # loads XCMSnExp object w. identified peaks called peaks_file
        load(paste0("./peaks_identified/", file)) 
        assign(name, peaks_file)
        if (i < length(files)) {
            concatenate_string <- paste0(c(concatenate_string, name, ","), collapse = "")
            }
        if (i == length(files)) {
            concatenate_string <- paste0(c(concatenate_string,  name, ")"), collapse = "")
            }
    }

    xset_peaksIdentified <- eval(parse(text = concatenate_string))
    save(xset_peaksIdentified, file = "./tmp/xset_peaksIdentified.Rdata")
    return(xset_peaksIdentified)
}
