
library(here)
library(GenomicRanges)

process_gff <- function(file) {
    ## There are some comments in the GFF file that contains useful
    ## information (e.g. seqlength).
    lines <- readLines(file)
    
    comments <- lines[grep("^#", lines)]
    comments <- comments[grep("^# Sequence Data", comments)] ## Also "Model Data"
    
    seq_lens <- stringr::str_remove(comments, "^# Sequence Data: ") %>%
        stringr::str_remove_all("seq...=") %>%
        stringr::str_remove_all('"') %>%
        strsplit(";") %>% sapply(function(lst) {
            seqnum <- as.integer(lst[[1]])
            seqlen <- as.integer(lst[[2]])
            seqhdr <- lst[[3]]
            names(seqlen) <- seqhdr
            seqlen
        })
    
    gff <- rtracklayer::readGFFAsGRanges(file)
    seqlengths(gff) <- seq_lens[names(seqlengths(gff))]
    gff
}

gff <- process_gff(here("sample_data/MAG06.gff"))
gff

seqlengths(gff)

li <- split(as.data.frame(gff), seqnames(gff), drop = TRUE)
li <- lapply(li, function(x) list(gene_track = (x)))
li <- lapply(li, function(x) {
    x$seqname <- unbox(as.character(unique(x$gene_track$seqnames)))
    x$seqlen <- unbox(as.integer(seqlengths(gff)[x$seqname]))
    x
})
li <- unname(li)
li[[1]]


dist <- here("srcweb/sample_data/")
if (!dir.exists(dist))
    dir.create(dist)

library(jsonlite)
for (el in li) {
    jsonlite::write_json(el, file.path(dist, paste0(el$seqname, ".json")), pretty = TRUE)
}



