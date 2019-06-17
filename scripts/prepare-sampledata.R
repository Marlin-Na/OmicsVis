
library(here)
library(GenomicRanges)

GFF <- here("sample_data/MAG06.gff")
EGGNOG <- here("sample_data/MAG06.faa.eggnog.txt")


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

process_eggnog <- function(file, gff) {
    eggnog <- readr::read_tsv(file, FALSE)
    id <- eggnog$X1
    eggnog_align_lst <- vector("list", length = length(gff))
    names(eggnog_align_lst) <- local({
        seqnames <- as.character(seqnames(gff))
        seqnames <- split(seqnames, seqnames)
        for (i in seq_along(seqnames)) {
            seqnames[[i]] <- paste0(seqnames[[i]], "_", seq(length(seqnames[[i]])))
        }
        unlist(seqnames)
    })
    eggnog_align_lst[id] <- split(eggnog, id)
    eggnog_align_lst[sapply(eggnog_align_lst, is.null)] <- list(data.frame())
    names(eggnog_align_lst) <- NULL
    eggnog_align_lst <- lapply(eggnog_align_lst, as.list)
    eggnog_align_lst
}

gff <- process_gff(GFF)
eggnog <- process_eggnog(EGGNOG, gff)

length(gff)
length(eggnog)

seqlengths(gff)

library(jsonlite)

df_gff <- as.data.frame(gff)
df_gff$eggnog <- eggnog
li <- split(df_gff, seqnames(gff), drop = TRUE)
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

for (el in li) {
    jsonlite::write_json(el, file.path(dist, paste0(el$seqname, ".json")), pretty = TRUE,
                         auto_unbox = TRUE)
}



