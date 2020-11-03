# input1 input2 output threads
args <- commandArgs(trailingOnly = TRUE)

library(dplyr)
library(tidyr)
library(data.table)
fixCollapsed <- function(df){
    colnames(df) <- c("key", "value")
    df <- df %>%
        mutate(key = strsplit(key, "; ")) %>%
        unnest(key)
    df <- df[, c(2, 1)]
    return(df)
}
fixDuplicated <- function(df){
    df <- df %>%
        group_by(key) %>%
        summarise(value = paste(value, collapse = "; "))
    values <- strsplit(df$value, "; ")
    values <- lapply(values, unique)
    values <- sapply(values, paste, collapse = "; ")
    df$value <- values
    return(df)
}
removeUnknown <- function(df){
    idx <- grepl("^-", df$key)
    df <- df[!idx,]
    return(df)
}
df <- fread("genbank2GO.txt", stringsAsFactors = FALSE,
    head = FALSE, nThread = parallel::detectCores())
df <- as.data.frame(df)
df %>%
    fixCollapsed() %>%
    fixDuplicated() %>%
    removeUnknown() %>%
    fwrite(file = "dictionary.txt", sep = "\t",
        nThread = parallel::detectCores())
df <- fread("refseq2GO.txt", stringsAsFactors = FALSE,
    head = FALSE, nThread = parallel::detectCores())
df <- as.data.frame(df)
df %>%
    fixCollapsed() %>%
    fixDuplicated() %>%
    removeUnknown() %>%
    fwrite(file = "dictionary.txt", sep = "\t",
        append = TRUE, nThread = parallel::detectCores())
