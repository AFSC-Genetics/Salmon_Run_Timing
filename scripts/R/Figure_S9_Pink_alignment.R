#!/usr/bin/env Rscript
# Simple PAF block plot: x = reference, y = haplotype contig, one segment per alignment.
# Genes are shaded along the x-axis. Uses only PAF columns 1,3,4,5,6,8,9 (no cs tag needed).

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# ---- parameters -------------------------------------------------------------
chr    <- "NC_060182.1"
rstart <- 72500000
rend   <- 73000000
flank  <- 25000          # extra reference shown each side; window is marked by dashed lines
flank  <- 25000          # extra reference shown each side; window is marked by dashed lines

paf_files <- c(hap1 = "Pink23E_4.bp.hap1.p_ctg.refpink.targetcontig.paf", hap2 = "Pink23E_4.bp.hap2.p_ctg.refpink.targetcontig.paf")
gff_file  <- "ncbi_dataset/data/pink_genome/genomic.gff"                            # full or subset GFF (.gz ok)

fstart <- rstart - flank
fend   <- rend + flank

# filter with awk so R only ever sees the rows it needs
open_cmd <- function(f) paste(if (grepl("\\.gz$", f)) "zcat" else "cat", shQuote(f))

# Reads only PAF columns 1,3,4,5,6,8,9 (the optional tag columns are ignored, so
# "parsing issues" warnings about extra columns are harmless). No shell commands.
read_blocks <- function(f, hap) {
  d <- read_tsv(f, col_names = FALSE, quote = "", progress = FALSE,
                col_types = cols_only(X1 = "c", X3 = "i", X4 = "i", X5 = "c",
                                      X6 = "c", X8 = "i", X9 = "i")) |>
    rename(qname = X1, qstart = X3, qend = X4, strand = X5,
           tname = X6, tstart = X8, tend = X9)
  hit <- d |>
    filter(tname == chr, tend >= fstart, tstart <= fend) |>   # any overlap with the plotted range
    select(-tname) |>
    mutate(hap = hap)
  if (nrow(hit) == 0)
    warning(hap, ": read ", nrow(d), " PAF lines but none overlap ", chr, ":", fstart, "-", fend,
            ". Target names seen: ", paste(head(unique(d$tname), 5), collapse = ", "))
  hit
}
blocks <- imap_dfr(paf_files, ~ read_blocks(.x, .y)) |>
  mutate(
    # segment from (tstart, q at tstart) to (tend, q at tend); minus strand runs qend -> qstart
    x0 = tstart, x1 = tend,
    y0 = if_else(strand == "+", qstart, qend),
    y1 = if_else(strand == "+", qend, qstart),
    # trim to the plotted range along the same line, so the y-axis stays sensible
    slope = (y1 - y0) / (x1 - x0),
    cx0 = pmax(x0, fstart), cx1 = pmin(x1, fend),
    cy0 = y0 + slope * (cx0 - x0),
    cy1 = y0 + slope * (cx1 - x0)
  )

# what was read: records per haplotype/contig, and how many extend past the plotted range
print(count(blocks, hap, qname, name = "n_alignments"))
message(sum(blocks$tstart < fstart | blocks$tend > fend), " of ", nrow(blocks),
        " alignments extend beyond the plotted range (kept, and trimmed at the edge for display)")

# ---- genes ------------------------------------------------------------------
gff_cmd <- sprintf("%s | awk -F'\\t' '$1==\"%s\" && ($3==\"gene\" || $3==\"pseudogene\") && $5>=%d && $4<=%d'",
                   open_cmd(gff_file), chr, as.integer(fstart), as.integer(fend))
genes <- read_tsv(pipe(gff_cmd),
                  col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"),
                  col_types = cols(.default = "c")) |>
  mutate(start = as.integer(start), end = as.integer(end),
         name = coalesce(str_extract(attr, "(?<=gene=)[^;]+"),
                         str_extract(attr, "(?<=Name=)[^;]+"),
                         str_extract(attr, "(?<=ID=)[^;]+")),
         xmin = pmax(start, fstart), xmax = pmin(end, fend))

# ---- plot -------------------------------------------------------------------
p <- ggplot() +
  geom_rect(data = genes, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = 0.4) +
  geom_segment(data = blocks, aes(x = cx0, xend = cx1, y = cy0, yend = cy1, colour = strand),
               linewidth = 1) +
  geom_text(data = genes, aes(x = (xmin + xmax) / 2, y = Inf, label = name),
            angle = 90, hjust = 1.05, size = 2.8, fontface = "italic", check_overlap = TRUE) +
  facet_grid(hap + qname ~ ., scales = "free_y", switch = "y") +
  scale_colour_manual(values = c("+" = "steelblue", "-" = "darkorange"), name = "Strand") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(fstart, fend), expand = FALSE) +
  labs(x = paste0("Reference ", chr, " position (bp)"),
       y = "Haplotype contig position (bp)") +
  theme_minimal(base_size = 11) +
  theme(strip.text.y.left = element_text(angle = 0), strip.placement = "outside",
        panel.grid.minor = element_blank())

n_panels <- nrow(distinct(blocks, hap, qname))
ggsave("paf_blocks.pdf", p, width = 9, height = 6)
ggsave("paf_blocks.png", p, width = 9, height = 6, dpi = 300)
