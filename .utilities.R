# .utilities.R
#
# Miscellaneous R code to suppport the Bioinformatics Introduction material
#
# Version: 1.2
# Date:    2017 01
# Author:  Boris Steipe
#
# V 1.2    added functions plotAAenrichment(), chunkSeq(),
#          dbAnnotateFeature(), dbFetchProteinData();
#          changed dbSanitizeSequence() to handle FASTA files.
# V 1.1    Simplify and remove course-specific code
# V 1.0    First code
#
# ToDo:
# Notes:
#
# ==============================================================================

objectInfo <- function(x) {
    # Function to combine various information items about R objects
    #
    # Input: an R object
    # Value: none - prints information as side-effect

    cat("object contents:")
    print(x, digits = 22)  # print value at maximal precision

    cat("\nstructure of object:\n")
    str(x)

    if (! is.list(x)) { # Don't use cat() if x is a list.
                        # cat() can't handle lists.
        cat("\nmode:   ", mode(x), "\n")
        cat("typeof: ", typeof(x), "\n")
        cat("class:  ", class(x), "\n")
    }

    # if the object has attributes, print them too
    if (! is.null(attributes(x))) {
        cat("\nattributes:\n")
        attributes(x)
    }
    # Done
}


plotAAenrichment <- function(s, ref) {

    # Purpose: Plot relative enrichment of amino acids.
    #
    # Parameters:
    #     s: A sequence in one-letter code
    #     ref: A vector of reference frequencies. If missing,
    #          aaindex[[459]] dataset from the seqinr package is used.
    #          If it has length 1, it is assumed to be a sequence, to
    #          be converted into a frequency table.
    # Value:
    #      The vector of sorted log2(fObs/fRef).
    #      Barplot as a sideeffect.
    #
    # Note: There are no sanity checks done on the input.
    #       There are also no checks that neither "ref" frequencies must not
    #       be zero for an amino acid in "obs".

    if (missing(ref)) {
        ref <- numeric()
        ref["A"] <- 7.9
        ref["R"] <- 4.9
        ref["N"] <- 4.0
        ref["D"] <- 5.5
        ref["C"] <- 1.9
        ref["Q"] <- 4.4
        ref["E"] <- 7.1
        ref["G"] <- 7.1
        ref["H"] <- 2.1
        ref["I"] <- 5.2
        ref["L"] <- 8.6
        ref["K"] <- 6.7
        ref["M"] <- 2.4
        ref["F"] <- 3.9
        ref["P"] <- 5.3
        ref["S"] <- 6.6
        ref["T"] <- 5.3
        ref["W"] <- 1.2
        ref["Y"] <- 3.1
        ref["V"] <- 6.8
    } else if (length(ref) == 1) {
        # assume this to be a sequence
        ref <- table(unlist(strsplit(toupper(ref), "")))
        ref <- 100 * (ref / sum(ref))  # Normalize
    }
    obs <- table(unlist(strsplit(toupper(s), "")))
    obs <- 100 * (obs / sum(obs))  # Normalize

    logR <- numeric()
    for (aa in names(obs)) {
        logR[aa] <- log(obs[aa] / ref[aa]) / log(2)
    }
    logR <- sort(logR, decreasing = TRUE)

    # plot
    label <- expression(paste(log[2],"( f(obs) / f(ref) )", sep = ""))

    chargePlus  <- "#282A8F"
    chargeMinus <- "#931523"
    hydrophilic <- "#57808F"
    hydrophobic <- "#6F7A79"
    plain       <- "#FDF7F7"

    # Preview colors
    # barplot(rep(1,5), col=c(chargePlus,
    #                         chargeMinus,
    #                         hydrophilic,
    #                         hydrophobic,
    #                         plain      ) )
    # Assign the colors to the different amino acid names
    barColors <- character()

    for (AA in names(logR)) {
        if      (grepl("[HKR]",     AA)) {barColors[AA] <- chargePlus }
        else if (grepl("[DE]",      AA)) {barColors[AA] <- chargeMinus}
        else if (grepl("[NQST]",    AA)) {barColors[AA] <- hydrophilic}
        else if (grepl("[FMILYVW]", AA)) {barColors[AA] <- hydrophobic}
        else                              barColors[AA] <- plain
    }

    barplot(logR,
            main = "AA composition enrichment",
            ylab = label,
            col = barColors,
            cex.names=0.8)
    abline(h=0)

    legend (x = 1,
            y = -1,
            legend = c("charged (+)",
                       "charged (-)",
                       "hydrophilic",
                       "hydrophobic",
                       "plain"),
            bty = "n",
            fill = c(chargePlus,
                     chargeMinus,
                     hydrophilic,
                     hydrophobic,
                     plain))
    # done

    invisible(logR)
}


biCode <- function(s) {

    # Purpose: make a 5 character code from a binomial name by concatenating
    #          the first three letter of the first word and the first
    #          two letters of the second word.
    #
    # Parameters:
    #     s: A scalar or vector of binomial species names.
    # Value:
    #      A scalar or vector of five-letter species codes.
    #
    # Note: there are no sanity checks done on the input.

    b <- character(length(s))
    for (i in 1:length(s)) {
        b[i] <- sprintf("%s%s",
                        toupper(unlist(substr(s[i], 1, 3))),
                        toupper(unlist(substr(strsplit(s[i], "\\s+")[[1]][2],
                                              1, 2))))
    }
    return(b)
}

dotPlot2 <- function(A, B,        # sequences
                     MDM = BLOSUM62,
                     f,           # filter
                     palette,     # a function that returns color values
                     xlab = "",
                     ylab = ""
) {
    # Purpose:
    #     Create a dotplot to measure sequence similarity between
    #     two amino acid sequences
    # Version:  1.0
    # Date:     2016-09
    # Author:   Boris Steipe
    #
    # Parameters:
    #     A, B: strings that contain no letters that are not found in MDM
    #     MDM: Mutation Data Matrix. Deafaults to BLOSUM62 which must be
    #          defined. Load package biostrings from BioConductor and load
    #          the matrix with data(BLOSUM62).
    #     f: filter matrix to weight an average around the neighborhood of
    #        an amino acid pair. Default to the identity matrix if missing.
    #        Average over a window of length f if length(f) is 1.
    #     palette: rainbow(), cm.colors(), or another function that returns
    #              a palette of color hexcodes. If missing, make our own
    #              palette.
    # Value:
    #     none. creates a dotplot.


    if (missing(f)) {
        f <- matrix(1) # default
    } else if (length(f) == 1) {
        if (! f %% 2) {stop("Sorry: f must be odd.")}
        w <- f
        f <- matrix(numeric(w * w), nrow = w)
        for (i in 1:w) { f[i, i] <- 1 }  # identity matrix
    }

    if (missing(palette)) {
        palette <- colorRampPalette(c("#000000",  #black
                                      "#111111",
                                      "#222222",
                                      "#332222",  # grey
                                      "#CC7755",
                                      "#EE8844",
                                      "#FF4400"), # red
                                    bias = 0.8)
    }
    A <- unlist(strsplit(A, ""))
    V <- unlist(strsplit(B, ""))
    lA <- length(A)
    lB <- length(B)

    m <- matrix(numeric(lA * lB), nrow = lA, ncol = lB)
    for (i in 1:lA) {
        for (j in 1:lB) {
            m[i, j] <- MDM[A[i], B[j]]
        }
    }
    m2 <- m
    wr <- floor((dim(f)[1] - 1) / 2)  # half-window size for rows
    wc <- floor((dim(f)[2] - 1) / 2)  # half-window size for columns

    for (i in (wr + 1):(lA - wr)) {
        for (j in (wc + 1):(lB - wc)) {
            # apply the filter to each value in m by weighting and summing
            # over its wr x wc neighborhood. Put the new value in m2
            m2[i, j] <- sum(f * m[(i-wr):(i+wr), (j-wc):(j+wc)])
        }
    }
    image(1:lA, 1:lB, m2,
          col = palette(24),
          ylim=c(lB,1), xlim=c(1,lA),
          xlab = xlab,
          ylab = ylab,
          axes = FALSE)
    box()

    # find good values for axis ticks and gridlines
    steps <- c(1, 2, 5, 10, 20, 50, 100, 200, 500,
               1000, 2000, 5000, 10000, 20000, 50000)
    gridStep <- sum(steps < max(lA, lB))

    # draw axes
    axis(1, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
    axis(2, at = c(1, seq(steps[gridStep - 3], lB, by=steps[gridStep - 3])))
    axis(3, at = c(1, seq(steps[gridStep - 3], lA, by=steps[gridStep - 3])))
    axis(4, at = c(1, seq(steps[gridStep - 2], lB, by=steps[gridStep - 3])))

    # draw grid with thin, transparent lines
    for (pos in seq(steps[gridStep - 2], lA, by = steps[gridStep - 2])) {
        abline(v=pos, col = "#FFFFFF44", lwd = 0.5)
    }
    for (pos in seq(steps[gridStep - 2], lB, by = steps[gridStep - 2])) {
        abline(h=pos, col = "#FFFFFF44", lwd = 0.5)
    }
}


pBar <- function(i, l, nCh = 50) {
    # Draw a progress bar in the console
    # i: the current iteration
    # l: the total number of iterations
    # nCh: width of the progress bar
    ticks <- round(seq(1, l-1, length.out = nCh))
    if (i < l) {
        if (any(i == ticks)) {
            p <- which(i == ticks)
            p1 <- paste(rep("#", p), collapse = "")
            p2 <- paste(rep("-", nCh - p), collapse = "")
            cat(sprintf("\r|%s%s|", p1, p2))
            flush.console()
        }
    }
    else { # done
        cat("\n")
    }
}



writeSeqSet <- function(seqSet,
                        file = stdout(),
                        format = "mfa",
                        blockSize = 50) {
    # Write a biostrings seqset to console or file in .mfa
    #  or .ali format.
    #
    # seqSet:    the MSA as MsaAAMultipleAlignment or AAStringSet object
    # file:      output filename. Default to console output.
    # format:    if format == "mfa", write a multi FASTA output
    #            if format == "ali", write a Claustal W output
    # blocksize: number of sequence characters per line.
    #
    MAXNAMEWID <- 15 # Maximum name width for "ali" format

    format <- str_trim(format)

    if (missing(seqSet)) {
        stop("Input object missing from arguments with no default.")
    }

    # Extract the raw data from the objects depending on
    # their respective class and put this
    # into a named vector of strings.
    if (class(seqSet)[1] == "MsaAAMultipleAlignment") {
        strings <- character(nrow(seqSet))
        for (i in 1:nrow(seqSet)) {
            strings[i] <- as.character(seqSet@unmasked[i])
            names(strings)[i] <- seqSet@unmasked@ranges@NAMES[i]
        }
    }
    else if (class(seqSet)[1] == "AAStringSet") {
        strings <- character(length(seqSet))
        for (i in 1:length(seqSet)) {
            strings[i] <- as.character(seqSet[i])
            names(strings)[i] <- seqSet@ranges@NAMES[i]
        }
    }
    else {
        stop(paste("Input object of class",
                   class(seqSet)[1],
                   "can't be handled by this function."))
    }

    out <- character()
    pattern <- sprintf(".{1,%d}", blockSize)

    if (format == "mfa") {
        for (i in 1:length(strings)) {
            # output FASTA header
            out <- c(out, sprintf(">%s", names(strings)[i]))
            # output sequence in blocks
            out <- c(out, unlist(str_match_all(strings[i], pattern)[[1]]))
        }
    } else if (format == "ali") {

        SEP <- paste(rep(" ", MAXNAMEWID), collapse="")
        out <- c(out, "CLUSTAL W formatted alignment", "")

        allBlocks <- character()
        for (i in 1:length(strings)) {
            # make labels for rownames
            label <- substr(paste(names(strings)[i],
                                  SEP,
                                  sep = ""), 1, MAXNAMEWID)
            label <- paste(label, " ", sep="")
            # chop strings into blocks
            blocks <- unlist(str_match_all(strings[i], pattern)[[1]])
            dim(blocks) <- c(1, length(blocks))
            rownames(blocks) <- label
            allBlocks <- rbind(allBlocks, blocks)
        }
        for (i in 1:ncol(allBlocks)) {
            for (j in 1:nrow(allBlocks)) {
                out <- c(out, sprintf("%s%s",
                                      rownames(allBlocks)[j],
                                      allBlocks[j, i]))
            }
            out <- c(out, "", "")
        }
    }
    writeLines(out, con = file)
}


maskSet <- function(set,
                    fGap = (2/3),
                    cGap="-",
                    verbose = TRUE) {
    # mask columns in "set" that contain more
    # then fGap fraction of cGap characters.

    if (class(set) != "MsaAAMultipleAlignment") {
        stop(paste("This function needs an object of class",
                   " MsaAAMultipleAlignment as input."))
    }

    lenAli <- set@unmasked@ranges@width[1]
    mat <- matrix(character(nrow(set) * lenAli),
                  ncol = lenAli)
    rownames(mat) <- set@unmasked@ranges@NAMES
    for (i in 1:nrow(set)) {
        seq <- as.character(set@unmasked[i])
        mat[i, ] <- unlist(strsplit(seq, ""))
    }
    colMask <- logical(lenAli)

    limit <- round(nrow(set) * fGap)
    for (i in 1:lenAli) {
        count <- table(mat[ , i])["-"]
        if (is.na(count)) { # No hyphen
            count <- 0
        }
        colMask[i] <- count <= limit
    }
    if (verbose) {
        cat(sprintf("Masking %4.2f %% of alignment columns.\n",
                    100 * (1 - (sum(colMask) / length(colMask)))))
    }
    mat <- mat[ , colMask]
    seqSet <- character()
    for (i in 1:nrow(mat)) {
        seqSet[i] <- paste(mat[i, ], collapse="")
    }
    names(seqSet) <- rownames(mat)
    return(AAStringSet(seqSet))
}

chunkSeq <- function(s, w = 50) {
    # Purpose: cat() sequence s in blocks of w characters.
    iChunks <- seq(1, nchar(s), by = w)
    chunks <- character(length(iChunks))
    cat("\n")
    for (i in 1:length(iChunks)) {
        chunks[i] <- substr(s, iChunks[i], iChunks[i] + w - 1)
        cat(chunks[i], "\n")
    }
    cat("\n")
    invisible(chunks)
}



# ====== DATABASE FUNCTIONS ============================================

dbInit <- function() {
    # Return an empty instance of the systems database
    # For the schema, see:
    # https://docs.google.com/presentation/d/1_nYWiwQc-9Z4anUAwOeVqWXgXIvM1Zl_AffMTM_No2Q

    db <- list()

    db$version <- "1.0"
    attributes(db$version)$tableCode <- "ver"

    db$system <- data.frame(
        ID = character(),
        name = character(),
        notes = character(),
        stringsAsFactors = FALSE)
    attributes(db$system$ID)$tableCode <- "sys"

    db$component <- data.frame(
        ID = character(),
        protein.ID = character(),
        system.ID = character(),
        status = character(),     # include / exclude / provisional
        notes = character(),
        stringsAsFactors = FALSE)
    attributes(db$component$ID)$tableCode <- "cmp"

    db$protein <- data.frame(
        ID = character(),
        name = character(),
        RefSeqID = character(),
        UniProtID = character(),
        taxonomy.ID = numeric(),
        sequence = character(),
        stringsAsFactors = FALSE)
    attributes(db$protein$ID)$tableCode <- "pro"

    db$taxonomy <- data.frame(
        ID = numeric(),
        species = character(),
        stringsAsFactors = FALSE)
    attributes(db$taxonomy$ID)$tableCode <- "tax"

    db$systemAnnotation <- data.frame(
        ID = character(),
        system.ID = character(),
        feature.ID = character(),
        stringsAsFactors = FALSE)
    attributes(db$systemAnnotation$ID)$tableCode <- "san"

    db$componentAnnotation <- data.frame(
        ID = character(),
        component.ID = character(),
        feature.ID = character(),
        stringsAsFactors = FALSE)
    attributes(db$componentAnnotation$ID)$tableCode <- "can"

    db$proteinAnnotation <- data.frame(
        ID = character(),
        protein.ID = character(),
        feature.ID = character(),
        start = numeric(),
        end = numeric(),
        stringsAsFactors = FALSE)
    attributes(db$proteinAnnotation$ID)$tableCode <- "pan"

    db$feature <- data.frame(
        ID = character(),
        name = character(),
        type.ID = character(),
        description = character(),
        sourceDB = character(),
        accession = character(),
        stringsAsFactors = FALSE)
    attributes(db$feature$ID)$tableCode <- "ftr"

    # type is one of:
    #     reference e.g. literature reference
    #     source (of information) e.g. literature reference
    #     3D domain e.g
    #     ptm       e.g. phosphorylation
    #     feature   e.g. low-complexity rtegion

    db$type <- data.frame(
        ID = character(),
        name = character(),
        description = character(),
        stringsAsFactors = FALSE)
    attributes(db$type$ID)$tableCode <- "typ"

    return(db)
}


dbSanitizeSequence <- function(s, strictAA = TRUE) {
    # Sanity check a protein sequence and convert to uppercase.
    #
    #        s:  an object that is or can be coerced to character
    # strictAA: modifies behavior regarding the presence of non one-letter
    #           code characters.
    #
    #    Value: a single string that contains only the 20 proteinogenic
    #           amino acids (if strictAA is TRUE) in uppercase, or
    #           other letters too (if strictAA is FALSE.) If s contains lines
    #           that are prefixed with a ">" (i.e. FASTA format), those
    #           header lines are removed.

    # Flatten any structure that s may have.
    s <- paste(unlist(s), collapse="")

    # It might be a FASTA format:
    if (substr(s, 1, 1) == ">") {
        # Split on linebreaks and remove
        # any lines beginning with ">"
        s <- unlist(strsplit(s, "\n"))
        s <- s[-(grep("^>", s))]
        s <- paste(s, collapse="")
    }

    # Remove all non-letters and convert to upper case.
    s <- toupper(gsub("[^a-zA-Z]", "", s))

    # Check for presence of non-amino acid letters.
    if (strictAA) {
        pattern <- "([BJOUXZ])"  # parentheses capture the match
        nonAA <- unlist(regmatches(s, regexec(pattern, s)))[1]
        if (! is.na(nonAA)) {
            stop(paste("Input contains letter outside of one-letter code: \"",
                       nonAA, "\".", sep=""))
        }
    }

    return(s)
}


dbConfirmUnique <- function(x) {
    # Confirm that x implements a unique selection on a table column.
    # x: a vector of logicals.
    # Value: x, if x contains exactly one TRUE element,
    #        stop() otherwise.

    if (any(!is.logical(x))) {
        stop("PANIC: Input fails is.logical() test.")
    } else if (sum(x) == 0) {
        stop("PANIC: No TRUE value in input.")
    } else  if (sum(x) > 1) {
        stop("PANIC: More than one TRUE value in input.")
    } else {
        return(x)
    }
}

dbAnnotateFeature <- function(db, pID, fID, Start, End) {
    # Purpose: Add an annotation to the proteinAnnotation table
    #    db: the target database
    #   pID: the protein ID
    #   fID: the feature ID
    # Start: start of the annotated range
    #   End: end of the annotated range
    # Value: the updated proteinAnnotation table

    panRow <- data.frame(ID = dbCreateID(db$proteinAnnotation$ID),
                         protein.ID = pID,
                         feature.ID = fID,
                         start = Start,
                         end = End,
                         stringsAsFactors = FALSE)
    return(rbind(db$proteinAnnotation, panRow))
}


dbCreateID <- function(x, ns = "my", tableCode) {

    # Purpose: Create a unique ID (primary key) in a databse table
    # Input:         x: a character vector of IDs structured as
    #                   <ns>_<tableCode>_<integer>
    #               ns: the "namespace" for the ID.
    #        tableCode: Optional. A short code to distinguish IDs in
    #                   different tables, for readability
    # Value:    a unique ID string <ns>_<tableCode>_<integer>, in which the
    #           integer part is one larger than the largest integer in that
    #           namespace in x.
    # Note:     tableCode is set as an attribute of the table ID
    #           column when the database first gets created. Thus it can be
    #           retrieved from x if that attribute is set. However we can't
    #           guarantee that it is always set because attributes don't
    #           survive rbind()'ing. Therefore we recover the tableCode from
    #           the first element, except if we are creating the firs ID
    #           in x.

    # tableCode ...
    if (missing(tableCode)) {
        if (length(x) == 0) {
            # make a tableCode for the first row in the table
            if (length(grep("tableCode", names(attributes(x)))) == 1) {
                tableCode <- attributes(x)$tableCode
            } else {
                # No attribute we can use, so create a default
                tableCode <- "xxx"
            }
        } else {
            # get tableCode from the first ID in x
            tableCode <- strsplit(x[1], "_")[[1]][2]
        }
    }

    # ID ...
    iIDs <- grep(sprintf("^%s_", ns), x)  # get indices of existing IDs for
                                           # requested namespace in x
    if (length(iIDs) == 0) { # no ID in this ns yet,
                              # return an ID with index "1"
        return(sprintf("%s_%s_%d", ns, tableCode, 1))
    } else {
        nums <- as.integer(gsub("^.+_", "", x[iIDs])) # remove prefix and
                                                      # coerce remaining
                                                      # strings to integer
        if (any(is.na(nums))) { # NA is present if coercion failed
            stop("PANIC: Input contains malformed ID(s).")
        } else {
            return(sprintf("%s_%s_%d",
                           ns,
                           tableCode,
                           as.integer(max(nums) + 1)))
        }
    }
}


dbGetFeatureSequence <- function(DB, name, feature) {
    # extracts a protein sequence from DB that matches the
    # protein name and the feature name, if these two are
    # unique, and has the start and end coordinates of the
    # feature annotation.

    proID <- DB$protein$ID[myDB$protein$name == name]
    if (length(proID) != 1) {stop("name does not match exactly one row.")}
    ftrID <- DB$feature$ID[myDB$feature$name == feature]
    if (length(ftrID) != 1) {stop("feature does not match exactly one row.")}
    fanID <- DB$proteinAnnotation$ID[
        DB$proteinAnnotation$protein.ID == proID &
            DB$proteinAnnotation$feature.ID == ftrID]
    if (length(fanID) != 1) {stop("annotation does not match exactly one row.")}

    return(substr(DB$protein$sequence[DB$protein$ID == proID],
                  DB$proteinAnnotation$start[DB$proteinAnnotation$ID == fanID],
                  DB$proteinAnnotation$end[DB$proteinAnnotation$ID == fanID]))
}

dbFetchProteinData <- function(rID, verbose = TRUE) {
    # Retrieves data for a given protein refseq ID from
    # NCBI and UniProt. NCBI data is retrived via the
    # eutils API.
    # UniProt IDs are retrieved from the
    # ID mapping service via a POST call.
    #
    # cf. http://www.uniprot.org/help/programmatic_access
    # rID: a valid RefSeqID
    #
    # Value: Returns a list of dataframe rows:
    #        db$version: db version for format
    #        db$taxonomy: one row of taxonomy data
    #        db$protein: one row of protein data
    #
    # Note! The db$protein$ID is not set but need to be created by the caller
    #       after the function returns.

    # load packages:
    if (!require(stringr, quietly=TRUE)) {
        install.packages("stringr")
        library(stringr)
    }

    if (!require(httr, quietly=TRUE)) {
        install.packages("httr")
        library(httr)
    }

    if (!require(XML, quietly=TRUE)) {
        install.packages("XML")
        library(XML)
    }

    if (!require(reutils, quietly=TRUE)) {
        install.packages("reutils")
        library(reutils)
    }

    # get UniProt ID for refseq
    response <- POST("http://www.uniprot.org/mapping/",
                     body = list(from = "P_REFSEQ_AC",
                                 to = "ACC",
                                 format = "tab",
                                 query = rID))
    if (response$status_code == 200) { # 200: oK
        uID <- str_trim(unlist(strsplit(httr::content(response), "\\s+"))[4])
        if (is.na(uID)) {
            warning(paste("UniProt ID mapping service returned NA.",
                          "Check your RefSeqID."))
        }
    } else {
        uID <- NA
        warning(paste("UniProt ID mapping not available:",
                      "server returned status",
                      response$status_code))
    }


    # Fetch NCBI data
    eHist <- epost(esearch(rID, "protein"), "protein")
    x <- content(esummary( eHist, version = "1.0"), "parsed")
    taxID <- x$TaxId
    binom <- efetch(taxID, "taxonomy")$xmlValue("/TaxaSet/Taxon/ScientificName")
    fasta <- content(efetch(eHist, retmode = "text", rettype = "fasta"))
    AA <- dbSanitizeSequence(fasta)

    # try using the word directly preceding the [<species name>] in the Title
    # as the name
    v <- unlist(strsplit(x$Title, "\\s+"))
    iSp <- grep("\\[", v)
    if (length(iSp) > 0 && iSp > 1) {
        nam <- v[iSp - 1]
    } else {
        nam <- NA
    }

    db <- list()
    db$version <- "1.0"
    db$protein <- data.frame(
        ID = NA,
        name = nam,
        RefSeqID = rID,
        UniProtID = uID,
        taxonomy.ID = taxID,
        sequence = AA,
        stringsAsFactors = FALSE)
    db$taxonomy <- data.frame(
        ID = taxID,
        species = binom,
        stringsAsFactors = FALSE)

    if (verbose) {
        cat(paste("Data retrieved:\n"))
        cat(paste("  name:         ", db$protein$name,         "\n"))
        cat(paste("  refSeqID:     ", db$protein$RefSeqID,     "\n"))
        cat(paste("  uniProtID:    ", db$protein$UniProtID,    "\n"))
        cat(paste("  taxID:        ", db$protein$taxonomy.ID,  "\n"))
        cat(paste("  organismName: ", db$taxonomy$species,     "\n"))
        len <- nchar(db$protein$sequence)
        cat(paste("  seq:           ",
                  substr(db$protein$sequence, 1, 10),
                  " ... ",
                  substr(db$protein$sequence, len-10, len),
                  " (", len, " AA)",
                  "\n\n", sep=""))
    }

    return(db)
}


# ====== SUPPORT FUNCTIONS =====================================================

waitTimer <- function(t, nIntervals = 50) {
    # pause and wait for t seconds and display a progress bar as
    # you are waiting
    t <- as.numeric(t)

    if (t < 0.1) {return(invisible())}

    increment <- t / nIntervals

    bar <- "----:----|"  # One module for the progress bar:
    bar <- rep(bar, ceiling(nIntervals / 10))  # repeat,
    bar <- unlist(strsplit(bar, "")) # split into single characters,
    bar <- bar[1:nIntervals]  # truncate,
    bar <- paste(bar, collapse="") # and collapse.

    cat(sprintf("\nWaiting: |%s\n         |", bar))
    for (i in 1:(nIntervals - 1)) {
        Sys.sleep(increment)
        cat("=")
    }
    Sys.sleep(increment)
    cat("|\n\n")

    return(invisible())
}


# ====== DATA ==================================================================


# 10 species of fungi for reference analysis.
# http://steipe.biochemistry.utoronto.ca/abc/index.php/Reference_species_for_fungi
REFspecies <- c("Aspergillus nidulans",
                "Bipolaris oryzae",
                "Coprinopsis cinerea",
                "Cryptococcus neoformans",
                "Neurospora crassa",
                "Puccinia graminis",
                "Saccharomyces cerevisiae",
                "Schizosaccharomyces pombe",
                "Ustilago maydis",
                "Wallemia mellicola"
)

# load the reference database refDB  (cf. create_refDB.R)
load("data/refDB.RData")




# [END]
