make_test_GRanges <- function() {
    new("GRanges",
        seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
        elementMetadata = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

test_GenomicRanges_shift <- function()
{
    ## empty, reversibility, recycling 'x'
    gr <- make_test_GRanges()
    checkIdentical(shift(GRanges(), 10), GRanges())
    checkIdentical(gr, shift(shift(gr, 10), -10))
    x <- 1:2
    checkIdentical(start(shift(gr[1:4], x)), start(gr[1:4]) + x)

    ## no seqlength or circularity
    checkIdentical(start(gr) + 10L, start(shift(gr, 10)))
    checkIdentical(width(gr), width(shift(gr, 10)))
    gr <- GRanges("chrA", IRanges(20, 30))
    checkIdentical(IRanges(8, 18), ranges(shift(gr, -12)))
    checkIdentical(IRanges(98, 108), ranges(shift(gr, 78)))

    ## seqlength and circularity combos
    gr <- GRanges("chr1", IRanges(5, width=6))
    isCircular(gr) <- TRUE 
    checkIdentical(start(shift(gr, -10)), -5L)

    seqlengths(gr) <- 20 
    isCircular(gr) <- NA
    warn <- FALSE 
    res <- withCallingHandlers({
        shift(gr, -10) 
    }, warning=function(w) {
        warn <<- TRUE 
        invokeRestart("muffleWarning")
    })
    checkTrue(warn == TRUE)
    checkIdentical(start(res), -5L)

    isCircular(gr) <- FALSE 
    warn <- FALSE 
    res <- withCallingHandlers({
        shift(gr, -10) 
    }, warning=function(w) {
        warn <<- TRUE 
        invokeRestart("muffleWarning")
    })
    checkTrue(warn == TRUE)
    checkIdentical(start(res), -5L)
}

test_GenomicRanges_trim <- function()
{
    checkIdentical(trim(GRanges()), GRanges())

    ## no seqlengths
    gr <- make_test_GRanges()
    checkIdentical(trim(gr), gr)

    ## seqlengths, isCircular NA and FALSE
    seqlengths(gr) <- c(10, NA, 20)
    spos <- suppressWarnings(shift(gr, 5))
    tend <- end(trim(spos))
    checkIdentical(tend, c(10L, rep(15L, 3), 10L, 10L, rep(15L, 4)))
    isCircular(gr)["chr1"] <- FALSE
    spos <- suppressWarnings(shift(gr, 5))
    tend <- end(trim(spos))
    checkIdentical(tend, c(10L, rep(15L, 3), 10L, 10L, rep(15L, 4)))

    ## seqlengths, isCircular TRUE 
    gr <- make_test_GRanges()
    seqlengths(gr) <- c(10, NA, 20)
    isCircular(gr)["chr1"] <- TRUE 
    spos <- suppressWarnings(shift(gr, 5))
    tend <- end(trim(spos))
    checkIdentical(tend, end(spos))
    spos <- suppressWarnings(shift(gr, 15))
    tend <- end(trim(spos))
    checkIdentical(tend, c(rep(25L, 6), rep(20L, 4)))
    isCircular(gr)["chr3"] <- TRUE 
    spos <- suppressWarnings(shift(gr, 15))
    tend <- end(trim(spos))
    checkIdentical(tend, end(spos))
}

test_GenomicRanges_flank <- function()
{
    gr <- make_test_GRanges()
    flanked <- flank(gr, 10)
    checkIdentical(rep(10L, length(gr)), width(flanked))
    checkIdentical(ifelse(as.vector(strand(gr) != "-"),
                          start(gr) - 10L, end(gr) + 1L), start(flanked))
    flanked <- flank(gr, 10, FALSE)
    checkIdentical(rep(10L, length(gr)), width(flanked))
    checkIdentical(ifelse(as.vector(strand(gr) != "-"),
                          end(gr) + 1L, start(gr) - 10L), start(flanked))
}

test_GenomicRanges_promoters <- function()
{
    checkTrue(length(promoters(GRanges())) == 0)

    ## upstream / downstream
    gr <- GRanges("chr1", IRanges(c(5, 10), width=1), "+")
    target <- GRanges("chr1", IRanges(c(5, 10), width=0), "+")
    current <- promoters(gr, 0, 0)
    checkIdentical(target, current)
    strand(gr) <- c("+", "-")
    target <- IRanges(c(3, 11), width=2)
    current <- ranges(promoters(gr, 2, 0))
    checkIdentical(target, current)
    target <- IRanges(c(5, 9), width=2)
    current <- ranges(promoters(gr, 0, 2))
    checkIdentical(target, current)

    gr <- GRanges("chr1", IRanges(0, width=6), "+")
    target <- GRanges("chr1", IRanges(-3, 2), "+")
    current <- promoters(gr, 3, 3)
    checkIdentical(target, current)
    checkTrue(validObject(current) == TRUE)
    gr <- GRanges("chr1", IRanges(rep(10, 3), width=6), c("+", "-", "*"))
    target <- GRanges("chr1", IRanges(c(7, 13, 7), c(12, 18, 12)),
        c("+", "-", "*"))
    current <- suppressWarnings(promoters(gr, 3, 3))
    checkIdentical(target, current)

    ## treat "*" as "+" 
    gr <- GRanges("chr1", IRanges(5, width=6), "+")
    target <- GRanges("chr1", IRanges(2, 7), "+")
    current <- promoters(gr, 3, 3)
    checkIdentical(target, current)
    strand(gr) <- "*"
    strand(target) <- "*"
    current <- suppressWarnings(promoters(gr, 3, 3))
    checkIdentical(target, current)

    ## metadata
    gr <- GRanges("chr1", IRanges(0, width=6), names="A", strand="+", score=99)
    current <- promoters(gr, 3, 3)
    checkIdentical(mcols(gr), mcols(current)) 
    checkIdentical(names(gr), names(current))
    checkIdentical(seqinfo(gr), seqinfo(current))
} 

test_GenomicRanges_resize <- function()
{
    gr <- make_test_GRanges()
    checkException(resize(gr, 10, fix = "middle"), silent = TRUE)
    checkException(resize(gr, 10, fix = rep("end", 3)), silent = TRUE)
    resized <- resize(gr, 10)
    checkIdentical(rep(10L, length(gr)), width(resized))
    checkIdentical(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 1L), start(resized))
    checkIdentical(ranges(resize(gr, 10, fix = "center")),
                   IRanges(rep(1:5, each=2), width = 10,
                           names = head(letters, 10)))
    checkIdentical(ranges(resize(gr, 10, fix = c("start", "end"))),
                   IRanges(c(1L, 1L, 3L, 1L, 5L, 1L, 7L, 1L, 1L, 10L),
                           width = 10, names = head(letters, 10)))
} 

test_GenomicRanges_restrict <- function()
{
    gr <-  make_test_GRanges()
    st <- structure(c(4,5), names = c("chr1", "chr2"))
    en <-  structure(c(8,9), names = c("chr2", "chr3"))
    res <- restrict(gr, start = st, end = en)
    checkIdentical(mcols(gr), mcols(res))
    checkIdentical(seqnames(gr), seqnames(res))
    checkIdentical(seqinfo(gr), seqinfo(res))
    target <- IRanges(start=c(4, 5, 5, 5, 5, 6, 7, 8, 9, 10),
                      end = c(10, 8, 8, 8, 10, 10, 9, 9, 9, 9),
                      names=letters[1:10])
    checkIdentical(ranges(res), target)
}

