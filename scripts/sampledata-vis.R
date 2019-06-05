
library(magrittr)
library(here)
library(TnT)

eggnog_annot <- readr::read_tsv(here("sample_data/MAG06.faa.eggnog.txt"), col_names = FALSE)

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

li.contig <- as.list(split(gff, seqnames(gff)))

view_contig <- function(n, additional_track = FALSE) {
    gr <- li.contig[[n]]
    label <- names(li.contig[n])
    ir <- ranges(gr)
    mcols(ir) <- mcols(gr)
    limit <- c(0, seqlengths(li.contig[[n]])[[unique(seqnames(li.contig[[n]]))]])
    
    reftrack <- TnT::BlockTrack(IRanges(start = limit[1], end = limit[2]), height = 20,
                                color = "grey", label = label)
    track <- TnT::BlockTrack(ir, NULL, color = "yellow", height = 20)
    
    if (additional_track) {
        ftrack <- TnT::FeatureTrack(ir, label = NULL, names = NULL)
        return(
        TnTBoard(
            list(reftrack, track, ftrack), allow.drag = FALSE, coord.range = limit,
            view.range = GRanges("UnKnown", IRanges(1, limit[2])) * 0.75
        )
        )
    }
    
    TnTBoard(
        list(reftrack, track), allow.drag = FALSE, coord.range = limit,
        view.range = GRanges("UnKnown", IRanges(1, limit[2])) * 0.75
    )
}
view_contig(3)
view_contig(3, T)
view_contig(5, T)
view_contig(7, T)
view_contig(9, T)
view_contig(10, T)
view_contig(14, T)
view_contig(22, T)
view_contig(102, T)



