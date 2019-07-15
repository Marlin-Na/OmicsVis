
library(here)
library(stringr)
library(dplyr)
library(GenomicRanges)

GFF <- here("sample_data/MAG06.gff")
EGGNOG <- here("sample_data/MAG06.faa.eggnog.txt")
EGGNOG_POS <- here("sample_data/MAG06.diamond.blast6.tsv")

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
    stopifnot(identical(str_remove(seqnames(gff), ".*_"), str_remove(gff$ID, "_.*")))
    ## Modify the ID so that it matchs with the eggnog
    gff$ID <- paste0(str_remove(seqnames(gff), "_.*"), "_", gff$ID)
    gff
}

gff <- process_gff(GFF)
df_gff <- tibble::as_tibble(as.data.frame(gff))
names(df_gff) <- paste0("gene", "_", names(df_gff))
eggnog <- readr::read_tsv(EGGNOG, FALSE)
names(eggnog) <- paste0("eggnog_", names(eggnog))


data <- dplyr::left_join(df_gff, eggnog, by = c("gene_ID" = "eggnog_X1"))
stopifnot(nrow(data) == nrow(df_gff))


#### The positional information from the diamond search used in eggnog

eggnog_pos <- readr::read_tsv(EGGNOG_POS, col_names = FALSE)
names(eggnog_pos) <- paste0("eggnog_pos_", names(eggnog_pos))
stopifnot(all(eggnog_pos$eggnog_pos_X1 %in% data$gene_ID))

## Calculate the start and end on the contig

eggnog_pos_offset <- data$gene_start[match(eggnog_pos$eggnog_pos_X1, data$gene_ID)]
eggnog_pos <- dplyr::mutate(
    eggnog_pos,
    eggnog_pos_start = eggnog_pos_offset + 3L * eggnog_pos_X7,
    eggnog_pos_end = eggnog_pos_offset + 3L * eggnog_pos_X8
)

data <- dplyr::left_join(data, eggnog_pos, by = c("gene_ID" = "eggnog_pos_X1"))

## gene* and eggnog* is a one-to-one mapping (though some gene does not have
## eggnog*, i.e. NA values)
## They are shown in the gene_track

## gene* and eggnog_pos* is one-to-many mapping and eggnog_pos*
## will be shown in a different track. 

gene_tracks <- dplyr::select(data, starts_with("gene"), starts_with("eggnog_X"))
gene_tracks <- dplyr::distinct(gene_tracks)
gene_tracks

## Generating random metrics for gene_track
set.seed(42L)
for (i in seq(30)) {
    ## Random values from 1 to 99
    gene_tracks[[sprintf("metric_%s", i)]] <- sample(1:99, nrow(gene_tracks), TRUE)
    ## 20% - 30% of the genes are expected to have a metric
    gene_tracks[[sprintf("metric_%s", i)]][
        sample(seq(nrow(gene_tracks)), round(nrow(gene_tracks)*3L/4L), TRUE)] <- NA
}

diamond_tracks <- dplyr::select(data, seqnames = "gene_seqnames",
                                "gene_ID", "gene_strand", starts_with("eggnog_pos"))
stopifnot(nrow(dplyr::distinct(diamond_tracks)) == nrow(diamond_tracks))
diamond_tracks

metas <- dplyr::select(data, seqname = "gene_seqnames")
metas <- dplyr::mutate(metas, seqlen = seqlengths(gff)[metas$seqname])
metas <- dplyr::distinct(metas)
metas

li.gene_tracks <- split(gene_tracks, gene_tracks$gene_seqnames)
li.metas <- split(metas, metas$seqname)
li.diamond_tracks <- split(diamond_tracks, diamond_tracks$seqnames)
stopifnot(identical(names(li.metas), names(li.gene_tracks)))
stopifnot(identical(names(li.metas), names(li.diamond_tracks)))

## Ad hoc way to remove the NA rows
for (i in seq_along(li.diamond_tracks)) {
    df <- li.diamond_tracks[[i]]
    df <- dplyr::filter(df, !is.na(eggnog_pos_X2))
    li.diamond_tracks[[i]] <- df
}

dist <- here("srcweb/sample_data/")
if (!dir.exists(dist))
    dir.create(dist)

for (i in seq_along(li.metas)) {
    ans <- list()
    ans$meta <- as.list(li.metas[[i]])
    ans$gene_track <- li.gene_tracks[[i]]
    ans$diamond_track <- li.diamond_tracks[[i]]
    dist_file <- file.path(dist, paste0(ans$meta$seqname, ".json"))
    print(dist_file)
    
    ## Open a binary connection so that it writes unix line ending
    dist_file <- file(dist_file, "wb")
    jsonlite::write_json(ans, dist_file, pretty = TRUE, auto_unbox = TRUE)
    close(dist_file)
}




#df_gff <- as.data.frame(gff)
#df_gff$eggnog <- eggnog
#li <- split(df_gff, seqnames(gff), drop = TRUE)
#li <- lapply(li, function(x) list(gene_track = (x)))
#li <- lapply(li, function(x) {
#    x$seqname <- unbox(as.character(unique(x$gene_track$seqnames)))
#    x$seqlen <- unbox(as.integer(seqlengths(gff)[x$seqname]))
#    x
#})
#li <- unname(li)
#li[[1]]
#
#
#
#for (el in li) {
#    jsonlite::write_json(el, file.path(dist, paste0(el$seqname, ".json")), pretty = TRUE,
#                         auto_unbox = TRUE)
#}


## Create index.json

index_table <- dplyr::select(
    metas, contig = seqname, length = seqlen) %>%
    left_join(gene_tracks[,c("gene_seqnames", "gene_ID")],
                by = c("contig" = "gene_seqnames")) %>%
    group_by(contig, length) %>%
    summarise(number_gene = length(unique(gene_ID)))

#dplyr::bind_rows(
#    lapply(li, function(x) {
#        contig <- as.character(x$seqname)
#        length <- as.integer(x$seqlen)
#        number_gene <- NROW(x$gene_track)
#        data.frame(contig = contig, length = length, number_gene = number_gene,
#                   stringsAsFactors = FALSE)
#    })
#)

index_table
jsonlite::write_json(index_table, file.path(dist, "index.json"),
                     pretty = TRUE, auto_unbox = TRUE)









# process_eggnog <- function(file, gff) {
#     eggnog <- readr::read_tsv(file, FALSE)
#     eggnog
#     #id <- eggnog$X1
#     #eggnog_align_lst <- vector("list", length = length(gff))
#     #names(eggnog_align_lst) <- local({
#     #    seqnames <- as.character(seqnames(gff))
#     #    seqnames <- split(seqnames, seqnames)
#     #    for (i in seq_along(seqnames)) {
#     #        seqnames[[i]] <- paste0(seqnames[[i]], "_", seq(length(seqnames[[i]])))
#     #    }
#     #    unlist(seqnames)
#     #})
#     #eggnog_align_lst[id] <- split(eggnog, id)
#     #eggnog_align_lst[sapply(eggnog_align_lst, is.null)] <- list(data.frame())
#     #names(eggnog_align_lst) <- NULL
#     #eggnog_align_lst <- lapply(eggnog_align_lst, as.list)
#     #eggnog_align_lst
# }
#eggnog <- process_eggnog(EGGNOG, gff)
