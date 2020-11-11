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
df <- fread(args[2], stringsAsFactors = FALSE,
    head = FALSE, nThread = as.integer(args[4]))
df <- as.data.frame(df)
df %>%
    fixCollapsed() %>%
    fixDuplicated() %>%
    removeUnknown() %>%
    fwrite(file = args[1], sep = "\t",
        nThread = as.integer(args[4]))
df <- fread(args[3], stringsAsFactors = FALSE,
    head = FALSE, nThread = as.integer(args[4]))
df <- as.data.frame(df)
df %>%
    fixCollapsed() %>%
    fixDuplicated() %>%
    removeUnknown() %>%
    fwrite(file = args[1], sep = "\t",
        append = TRUE, nThread = as.integer(args[4]))
