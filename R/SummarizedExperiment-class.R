## 
## Helper classes: Assays, ShallowData, ShallowSimpleListAssays

setClass("Assays")

.ShallowData <- setRefClass("ShallowData",
    fields = list( data = "ANY" ))

setMethod("clone", "ShallowData",  # not exported
    function(x, ...)
{
    args <- list(...)
    x <- x$copy()
    for (nm in names(args))
        x$field(nm, args[[nm]])
    x
})

.ShallowSimpleListAssays0 <- setRefClass("ShallowSimpleListAssays",
    fields = list( data = "SimpleList" ),
    contains = c("ShallowData", "Assays"))

.ShallowSimpleListAssays <-
    function(..., data = SimpleList())
    ## this two-step constructor is much faster, 16 May, 2013
{
    xx <- .ShallowSimpleListAssays0(...)
    xx$data <- data
    xx
}

## 
## SummarizedExperiment

setClass("SummarizedExperiment",
    representation(
        exptData="SimpleList",                # overall description
        rowData="GenomicRangesORGRangesList", # rows and their annotations
        colData="DataFrame",                  # columns and their annotations
        assays="Assays"),                     # Data -- e.g., list of matricies
    prototype(
        rowData=GRanges(),
        assays=.ShallowSimpleListAssays(data=SimpleList())))

## validity

.valid.SummarizedExperiment.assays_current <- function(x)
{
    if (!is(slot(x, "assays"), "Assays"))
        return("'assays' is out-of-date; use updateObject()")
    NULL
}

.valid.SummarizedExperiment.assays_class <- function(x)
{
    ok <- sapply(assays(x, withDimnames=FALSE), function(cl) {
        (!is.null(dim(cl))) && (length(dim(cl)) >= 2L)
    })
    if (!all(ok))
        return("'assays' must be matrix-like with 2 (or more?) dimensions")
    NULL
}

.valid.SummarizedExperiment.rowData_dims <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), nrow) ==
             length(rowData(x))))
        return("'rowData' length differs from 'assays' nrow")
    NULL
}

.valid.SummarizedExperiment.colData_dims <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), ncol) ==
             nrow(colData(x))))
        return("'colData' nrow differs from 'assays' ncol")
    NULL
}

.valid.SummarizedExperiment.assays_dims <- function(x)
{
    c(.valid.SummarizedExperiment.rowData_dims(x),
      .valid.SummarizedExperiment.colData_dims(x))
}

.valid.SummarizedExperiment <- function(x)
{
    c(.valid.SummarizedExperiment.assays_current(x),
      msg <- .valid.SummarizedExperiment.assays_class(x),
      if (is.null(msg)) {
          .valid.SummarizedExperiment.assays_dims(x)
      } else NULL)
}

setValidity2("SummarizedExperiment", .valid.SummarizedExperiment)

setGeneric("SummarizedExperiment",
    function(assays, ...) standardGeneric("SummarizedExperiment"))

.GRangesList_assays <-
    function(assays)
{
    m <- assays[[1]]
    n <- nrow(m)
    names <- rownames(m)
    relist(GRanges(), PartitioningByEnd(integer(n), names=names))
}

setMethod(SummarizedExperiment, "SimpleList",
   function(assays, rowData=GRangesList(), colData=DataFrame(),
             exptData=SimpleList(), ...,
             verbose=FALSE)
{
    if (missing(rowData) && 0L != length(assays))
        rowData <- .GRangesList_assays(assays)
    if (missing(colData) && 0L != length(assays)) {
        nms <- colnames(assays[[1]])
        if (is.null(nms) && 0L != ncol(assays[[1]]))
            stop("'SummarizedExperiment' assay colnames must not be NULL")
        colData <- DataFrame(row.names=nms)
    }

    FUN <- function(x)
        ## dimnames as NULL or list(NULL, NULL)
        all(sapply(dimnames(x), is.null))
    if (!all(sapply(assays, FUN)))
        assays <- endoapply(assays, unname)
    if (!is(assays, "Assays"))
        assays <- .ShallowSimpleListAssays(data=assays)

    new("SummarizedExperiment", exptData=exptData, rowData=rowData,
        colData=colData, assays=assays, ...)
})

setMethod(SummarizedExperiment, "missing",
    function(assays, ...)
{
    SummarizedExperiment(SimpleList(), ...)
})

setMethod(SummarizedExperiment, "list",
    function(assays, ...)
{
    SummarizedExperiment(do.call(SimpleList, assays), ...)
})

setMethod(SummarizedExperiment, "matrix",
    function(assays, ...)
{
    SummarizedExperiment(SimpleList(assays), ...)
})

## update / clone

setMethod("clone", "SummarizedExperiment",  # not exported
    function(x, ...)
{
    ## IRanges:::extraArgsAsList would prevent using clone on
    ## subclasses
    args <- list(...)
    firstTime <- TRUE
    for (nm in names(args)) {
        s <- slot(x, nm)
        v <- args[[nm]]
        if (is(s, "ShallowData"))
            v <- clone(s, data=v)
        if (firstTime) {
            slot(x, nm, FALSE) <- v
            firstTime <- FALSE
        } else {
            `slot<-`(x, nm, FALSE, v)
        }
    }
    x
})

setGeneric("value",  # not exported
    function(x, name, ...) standardGeneric("value"),
    signature = "x")

setMethod("value", "SummarizedExperiment", # not exported
    function(x, name, ...)
{
    s <- slot(x, name)
    if (is(s, "ShallowData"))
        s <- s$data
    s
})

## getter / setter generics

setGeneric("exptData", function(x, ...) standardGeneric("exptData"))

setGeneric("exptData<-",
    function(x, ..., value) standardGeneric("exptData<-"))

## rowData, colData seem too vague, but from eSet derived classes wanted to
## call the rows / cols something different from 'features' or 'samples', so
## might as well avoid the issue
setGeneric("rowData", function(x, ...) standardGeneric("rowData"))

setGeneric("rowData<-",
    function(x, ..., value) standardGeneric("rowData<-"))

setGeneric("colData", function(x, ...) standardGeneric("colData"))

setGeneric("colData<-",
    function(x, ..., value) standardGeneric("colData<-"))

setGeneric("assays",
    function(x, ..., withDimnames=TRUE) standardGeneric("assays"),
    signature="x")

setGeneric("assays<-",
    function(x, ..., value) standardGeneric("assays<-"))

setGeneric("assay", function(x, i, ...) standardGeneric("assay"))

setGeneric("assay<-",
    function(x, i, ..., value) standardGeneric("assay<-"))

## Simple 'getters' / 'setters'

setMethod(exptData, "SummarizedExperiment",
    function(x, ...) value(x, "exptData"))

setReplaceMethod("exptData", c("SummarizedExperiment", "SimpleList"),
    function(x, ..., value)
{
    clone(x, ..., exptData=value)
})

setReplaceMethod("exptData", c("SummarizedExperiment", "list"),
    function(x, ..., value)
{
    clone(x, ..., exptData=SimpleList(value))
})

setMethod(rowData, "SummarizedExperiment",
    function(x, ...) value(x, "rowData"))

.SummarizedExperiment.rowData.replace <-
    function(x, ..., value)
{
    x <- clone(x, ..., rowData=value)
    msg <- .valid.SummarizedExperiment.rowData_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("rowData", c("SummarizedExperiment", "GenomicRanges"),
    .SummarizedExperiment.rowData.replace)

setReplaceMethod("rowData", c("SummarizedExperiment", "GRangesList"),
    .SummarizedExperiment.rowData.replace)

setMethod(colData, "SummarizedExperiment",
    function(x, ...) value(x, "colData"))

setReplaceMethod("colData", c("SummarizedExperiment", "DataFrame"),
    function(x, ..., value)
{
    x <- clone(x, ..., colData=value)
    msg <- .valid.SummarizedExperiment.colData_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
})

setMethod(assays, "SummarizedExperiment",
    function(x, ..., withDimnames=TRUE)
{
    if (withDimnames)
        endoapply(value(x, "assays"), "dimnames<-", dimnames(x))
    else
        value(x, "assays")
})

.SummarizedExperiment.assays.replace <-
    function(x, ..., value)
{
    x <- clone(x, ..., assays=value)
    msg <- .valid.SummarizedExperiment(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("assays", c("SummarizedExperiment", "SimpleList"),
    .SummarizedExperiment.assays.replace)

setReplaceMethod("assays", c("SummarizedExperiment", "list"),
    function(x, ..., value)
{
    value <- lapply(x, "dimnames<-", NULL)
    .SummarizedExperiment.assays.replace(x, ...,
                                         value=SimpleList(value))
})

## convenience for common use case 
setMethod(assay, c("SummarizedExperiment", "missing"),
    function(x, i, ...)
{
    assays <- assays(x, ...)
    if (0L == length(assays))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(assays(<", class(x), ">)) is 0'")
    assays[[1]]
})

setMethod(assay, c("SummarizedExperiment", "numeric"),
    function(x, i, ...) 
{
    tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

setMethod(assay, c("SummarizedExperiment", "character"),
    function(x, i = names(x)[1], ...) 
{
    msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                  "invalid subscript 'i'")
    res <- tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
    if (is.null(res))
        stop(msg, "\n'i' not in names(assays(<", class(x), ">))")
    res
})

setReplaceMethod("assay", c("SummarizedExperiment", "missing", "matrix"),
    function(x, i, ..., value)
{
    if (0L == length(assays(x)))
        stop("'assay(<", class(x), ">) <- value' ", "length(assays(<",
             class(x), ">)) is 0")
    assays(x)[[1]] <- value
    x
})

setReplaceMethod("assay",
    c("SummarizedExperiment", "numeric", "matrix"),
    function(x, i = 1, ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

setReplaceMethod("assay",
    c("SummarizedExperiment", "character", "matrix"),
    function(x, i = names(x)[1], ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

## cannonical location for dim, dimnames; dimnames should be checked
## for consistency (if non-null) and stripped from assays on
## construction, or added from assays if row/col names are NULL in
## <SummarizedExperiment> but not assays. dimnames need to be added on
## to assays when assays() invoked
setMethod(dim, "SummarizedExperiment",
    function(x) 
{
    c(length(rowData(x)), nrow(colData(x)))
})

setMethod(dimnames, "SummarizedExperiment",
    function(x)
{
    list(names(rowData(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("SummarizedExperiment", "list"),
    function(x, value) 
{
    rowData <- rowData(x)
    names(rowData) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    clone(x, rowData=rowData, colData=colData)
})
 
setReplaceMethod("dimnames", c("SummarizedExperiment", "NULL"),
    function(x, value)
{
    dimnames(x) <- list(NULL, NULL)
    x
})

## Subset -- array-like

.SummarizedExperiment.charbound <-
    function(idx, txt, fmt)
{
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
        msg <- paste(IRanges:::selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

.SummarizedExperiment.assays.subset <-
    function(x, i, j, ...)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x) {
            if (length(dim(x)) == 2L)
                x[i, j, drop=FALSE]
            else
                x[i, j, , drop=FALSE]
        }
    } else if (!missing(i)) {
        fun <- function(x) {
            if (length(dim(x)) == 2L)
                x[i, , drop=FALSE]
            else
                x[i, , , drop=FALSE]
        }
    } else if (!missing(j)) {
        fun <- function(x) {
            if (length(dim(x)) == 2L)
                x[, j, drop=FALSE]
            else
                x[, j, , drop=FALSE]
        }
    }
    endoapply(assays(x, withDimnames=FALSE), fun)
}

setMethod("[", c("SummarizedExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    if (missing(i) && missing(j))
        return(x)

    if (!missing(i) && is.character(i)) {
        fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
        i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
    }
    if (!missing(j) && is.character(j)) {
        fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
        j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
    }

    if (!missing(i) && !missing(j)) {
        ii <- as.vector(i)
        jj <- as.vector(j)
        x <- clone(x, ..., rowData=rowData(x)[i],
            colData=colData(x)[j, , drop=FALSE],
            assays=.SummarizedExperiment.assays.subset(x, ii, jj))
    } else if (missing(i)) {
        jj <- as.vector(j)
        x <- clone(x, ..., colData=colData(x)[j, , drop=FALSE],
            assays=.SummarizedExperiment.assays.subset(x, j=jj))
    } else {                            # missing(j)
        ii <- as.vector(i)
        x <- clone(x, ..., rowData=rowData(x)[i],
            assays=.SummarizedExperiment.assays.subset(x, ii))
    }
    x
})

.SummarizedExperiment.assays.subsetgets <-
    function(x, i, j, ..., value)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x, value) {
            if (length(dim(x)) == 2L)
                x[i, j] <- value
            else
                x[i, j, ] <- value
            x
        }
    } else if (!missing(i)) {
        fun <- function(x, value) {
            if (length(dim(x)) == 2L)
                x[i, ] <- value
            else
                x[i, , ] <- value
            x
        }
    } else if (!missing(j)) {
        fun <- function(x, value) {
            if (length(dim(x)) == 2L)
                x[, j] <- value
            else
                x[, j, ] <- value
            x
        }
    }
    a <- assays(x, withDimnames=FALSE)
    v <- assays(value, withDimnames=FALSE)
    mendoapply(fun, x=a, value=v)
}

setReplaceMethod("[",
    c("SummarizedExperiment", "ANY", "ANY", "SummarizedExperiment"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        return(value)
 
    if (!missing(i) && is.character(i)) {
        fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
        i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
    }

    if (!missing(j) && is.character(j)) {
        fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
        j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
    }

    if (!missing(i) && !missing(j)) {
        ii <- as.vector(i)
        jj <- as.vector(j)
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            rowData=local({
                r <- rowData(x)
                r[i] <- rowData(value)
                names(r)[ii] <- names(rowData(value))
                r
            }), colData=local({
                c <- colData(x)
                c[j,] <- colData(value)
                rownames(c)[jj] <- rownames(colData(value))
                c
            }), assays=.SummarizedExperiment.assays.subsetgets(x, ii, jj,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.assays_dims(x)
    } else if (missing(i)) {
        jj <- as.vector(j)
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            colData=local({
                c <- colData(x)
                c[j,] <- colData(value)
                rownames(c)[jj] <- rownames(colData(value))
                c
            }), assays=.SummarizedExperiment.assays.subsetgets(x, j=jj,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.colData_dims(x)
    } else {                            # missing(j)
        ii <- as.vector(i)
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            rowData=local({
                r <- rowData(x)
                r[i] <- rowData(value)
                names(r)[ii] <- names(rowData(value))
                r
            }), assays=.SummarizedExperiment.assays.subsetgets(x, ii,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.rowData_dims(x)
    }
    if (!is.null(msg))
        stop(msg)
    x
})

## rbind, cbind

## Appropriate for objects with different ranges and same samples.
setMethod("rbind", "SummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .rbind.SummarizedExperiment(args)
})

.rbind.SummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, colnames)))
            stop("'...' objects must have the same colnames")
    if (!.compare(lapply(args, ncol)))
            stop("'...' objects must have the same number of samples")
    rowData <- do.call(c, lapply(args, 
        function(i) slot(i, "rowData")))
    cdlst <- lapply(args, colData) 
    colData <- .cbind.DataFrame(args, colData, "colData") 
    assays <- .bind.arrays(args, rbind, "assays")
    exptData <- do.call(c, lapply(args, exptData))

    SummarizedExperiment(assays=assays, rowData=rowData, 
                         colData=colData, exptData=exptData) 
}

## Appropriate for objects with same ranges and different samples.
setMethod("cbind", "SummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .cbind.SummarizedExperiment(args)
})

.cbind.SummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, rowData), TRUE))
        stop("'...' object ranges (rows) are not compatible")
    rowData <- rowData(args[[1]])
    mcols(rowData) <- .cbind.DataFrame(args, mcols, "mcols")
    colData <- do.call(rbind, lapply(args, colData))
    assays <- .bind.arrays(args, cbind, "assays")
    exptData <- do.call(c, lapply(args, exptData))

    SummarizedExperiment(assays=assays, rowData=rowData, 
                         colData=colData, exptData=exptData) 
}

.compare <- function(x, GenomicRanges=FALSE) 
{
    x1 <- x[[1]]
    if (GenomicRanges) {
        if (is(x1, "GRangesList")) {
            x <- lapply(x, unlist) 
            x1 <- x[[1]]
        }
        for (i in seq_along(x[-1])) {
            if (length(x1) != length(x[[i]]))
                return(FALSE)
            ok <- x1 == x[[i]]
            if (!all(ok))
                return(FALSE)
        }
        return(TRUE) 
    } else {
        all(sapply(x[-1], 
            function(xelt) all(identical(xelt, x[[1]]))))
    }
}

.cbind.DataFrame <- function(args, accessor, accessorName)
{
    lst <- lapply(args, accessor)
    if (!.compare(lst)) {
        nms <- lapply(lst, names)
        nmsv <- unlist(nms, use.names=FALSE)
        names(nmsv) <- rep(seq_along(nms), elementLengths(nms))
        dups <- duplicated(nmsv)
        ## no duplicates 
        if (!any(dups))
            return(do.call(cbind, lst))
        ## confirm duplicates are the same 
        lapply(nmsv[duplicated(nmsv)], function(d) {
            if (!.compare(lapply(lst, "[", d)))
                stop("column(s) '", unname(d),
                     "' in ", sQuote(accessorName),
                     " are duplicated and the data do not match")})
        ## remove duplicates
        do.call(cbind, lst)[,!dups]
    } else {
        lst[[1]]
    }
}

.bind.arrays <- function(args, bind, accessor, accessorName)
{
    lst <- lapply(args, accessor)
    var <- lapply(lst,  names)
    if (is.null(uvar <- unique(unlist(var))))
        return(SimpleList())
    if (!.compare(var))
        stop("elements in ", sQuote(accessorName), 
             " must have the same names")

    ## extract variable
    l1 <- lapply(uvar, function(v) {
        l2 <- lapply(lst, "[[", v)
        dim <- .assay.dimension(l2, bind)
        if (!is.na(dim[3])) {
            ## assay: combine each dimension of each element
            l3 <- lapply(1:dim[3], 
                      function(i) do.call(bind, lapply(l2, "[", ,,i)))
            array(do.call(c, l3), dim=dim)
        } else {
            ## matrix: combine each element
            do.call(bind, l2)
        }
    })
    names(l1) <- uvar
    SimpleList(l1)
}

.assay.dimension <- function(lst, bind)
{
    dim <- lapply(lst, dim)
    if (!.compare(lapply(dim, "[", 3)))
        stop("elements in assays must have the same dimension")
    if (identical(bind, cbind))
        c(dim[[1]][1], do.call(sum, lapply(dim, "[", 2)), dim[[1]][3])
    else
        c(do.call(sum, lapply(dim, "[", 1)), dim[[1]][2], dim[[1]][3])
}

## $, $<-, [[, [[<- for colData access

setMethod("[[", c("SummarizedExperiment", "ANY", "missing"),
    function(x, i, j, ...)
{
    colData(x)[[i, ...]]
})

setReplaceMethod("[[",
    c("SummarizedExperiment", "ANY", "missing", "ANY"),
    function(x, i, j, ..., value)
{
    colData(x)[[i, ...]] <- value
    x
})

setMethod("$", "SummarizedExperiment",
    function(x, name)
{
    colData(x)[[name]]
})

setReplaceMethod("$", c("SummarizedExperiment", "ANY"),
    function(x, name, value)
{
    colData(x)[[name]] <- value
    x
})

## compatibility

setMethod(updateObject, "SummarizedExperiment",
    function(object, ..., verbose=FALSE)
{
    s <- slot(object, "assays")
    if (is(s, "SimpleList"))
        slot(object, "assays") <- .ShallowSimpleListAssays(data=s)
    object
})

## show

setMethod(show, "SummarizedExperiment",
    function(object)
{
    selectSome <- IRanges:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(IRanges:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    expt <- names(exptData(object))
    if (is.null(expt))
        expt <- character(length(exptData(object)))
    scat("exptData(%d): %s\n", expt)
    nms <- names(assays(object, withDimnames=FALSE))
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    scat("assays(%d): %s\n", nms)
    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]]) scat("rownames(%d): %s\n", dimnames[[1]])
    else scat("rownames: NULL\n")
    scat("rowData metadata column names(%d): %s\n",
         names(mcols(rowData(object))))
    if (dlen[[2]]) scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")
    scat("colData names(%d): %s\n", names(colData(object)))
})
