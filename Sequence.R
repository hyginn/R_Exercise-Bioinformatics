# Sequence.R
#
# Purpose:  Introduction to Bioinformatics with R - Sequence
#
# Version: 1.0
#
# Date:    2017  01  03
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First code
#
# TODO:
#
#
# ==============================================================================

# ==============================================================================
#      PART ONE: ADDING DATA
# ==============================================================================

# To add data for SPIPU, let's first make a copy of refDB, call it myDB
myDB <- refDB

# Add the MBP1_SPIPU sequence:
iRow <- nrow(myDB$protein) + 1

# Confirm that all of the data items we enter below are the same that you got in
# your NCBI and Uniprot searches:
#
# protein table:
( myDB$protein[iRow, "ID"] <- dbCreateID(myDB$protein$ID) )
myDB$protein[iRow, "name"] <- "MBP1_SPIPU"
myDB$protein[iRow, "RefSeqID"] <- "XP_016612325"
myDB$protein[iRow, "UniProtID"] <- "A0A0L0HSF4"
myDB$protein[iRow, "taxonomy.ID"] <- 645134
myDB$protein[iRow, "sequence"] <- dbSanitizeSequence(
" 1 metppsvpvs vpangssggg ptiysaiysg vpvyemmcrd tpvmrrmada ymnatqilkv
 61 agftkpkrtk ilerevlqge hekiqggygk yqgtwiplse srmlasrynv dtelrplfdf
121 dpqtgsvapk pptlrtphvd lagdglasps spsppphhsl dyfspvdgps hllsefnspd
181 gsitpqgtvs pppgineeap qptwrdnlnl ssmtreassq pldqphlnsq vqsilnsfgm
241 qhllqhpstd gaisntesss nantvnemet ptttnsqetn rqtsssqqsg ptpssrptqs
301 qatsieaeyy qklqqgllld estaqnyvah iqqilrtqdm qqfqqaqqys sylqagsypa
361 vqsayqqlam lnllslaaqh aytypngyps ldswfkqtdg taskslsglg qpitgqrsqv
421 ssspeqskvp sndrrqsgve edhdapelds rhtsahqlrq revltnlyly sdsqlhkiid
481 slrnptnprk ldvnlvlegk rgstllhfaa algreklvka lisnganpai aanggetplm
541 ravyhvgcyd lqtfsqiadw lrdsvyltdt srrtvlhhia kaskkrdlwa aagyyvrcia
601 qviqeeysea knlipqvvna qdaygntalh yackyasrsv aevllklgat heianhaget
661 pvdlaryepk llallfreie rddvaawddf dedgdilmrd fdeatdtddd eferspspnp
721 versmheqtm ntldkfygek venilsnltp svtnavkagv dqlrnnlmsr yreqearvei
781 tqleleqake qlanakkels elqekrkrmd elvelanele dkleqhatrl irktkgkrav
841 epsissdedk edddthdpdm vgeeatqsgp saeaesaarp leddpskpsr dpsdsarmnn
901 eaipsqstdp essisqtstt krkrsasltd dedtvsdaaa hssdtshtsk lrklmpsdrt
961 rtrsgrqaev rlkaeladlq qeverseqnr lkmfheieni rahreekdrk fkriiaacvr
1021 mpmdkideli eplvhele")

# taxonomy table:
iRow <- nrow(myDB$taxonomy) + 1
myDB$taxonomy[iRow, "ID"] <- 645134
myDB$taxonomy[iRow, "species"] <- "Spizellomyces punctatus DAOM BR117"

# Confirm: Is our new protein named according to the naming pattern of the
# others? It should be. And does the taxonomy table contain the binomial name?
# And what about its features? Let's compute sequence lengths on the fly (with
# the function nchar() ) and open a table of selected information with the table
# viewer function View().

View(cbind(myDB$protein[ , c("ID", "name", "RefSeqID")],
           length = nchar(myDB$protein$sequence)))

# Where does our new protein's length fall relative to the reference proteins?
# About the same? Much shorter? Much longer?

# Is that the right sequence?
sel <- myDB$protein$ID == "my_pro_1"
myDB$protein$sequence[sel]

# Is that the right taxonomy ID and binomial name for YFO?
sel <- biCode(myDB$taxonomy$species) == "SPIPU"
myDB$taxonomy[sel, ]

# If there are any errors or inconsistencies, do not continue, fix the problem
# first.

# === Saving and loading data ==================================================
#

# Once you have confirmed that myDB has been successfully created and the
# MBP1_SPIPU protein added, and the database is in a good state, you should make
# a backup copy. There are many ways to save data to a file on disk and read it
# back in. One of the most convenient is the function pair save() and load().
#
# save() saves an R object to a file. Its signature is

#    save(<object-name or names>, file = <file-name>)

# The object, or objects that are saved like this are identically recreated when
# load() is executed. The object name is not changed - you don't assign the
# result of load(). You use load() for its "side-effect": re-creating the saved
# object, not for using the return-value. Therefore the signature for load() is
# simply

#     load(<file-name>)

# All R objects in <file-name> are created by load(), if they don't yet exist.
# If they already exist, the will be overwritten. The only effect you see is
# that the object appears in the Environment pane of RStudio (if it wasn't there
# already), or you may notice that its content has changed if it did exist.

# == TASK: ==
# Save myDB so you can recreate it easily when you next open RStudio.

save(myDB, file = "myDB.01.RData")  # Note that I give this a version number!

# Let's confirm:
rm(myDB)
nrow(myDB$protein)  # must give an error
load("myDB.01.RData")
nrow(myDB$protein)  # must be 11. If not, don't continue etc.



# ==============================================================================
#      PART TWO: REGULAR EXPRESSIONS
# ==============================================================================

# We have encountered regular expressions briefly in the Data.R script, where we
# used them in our dbSanitizeSequence() function. Let's study them a bit more
# deeply now.

# === STRING MATCHING

# Not all pattern matching uses (and needs) the capability of regular
# expressions to specify variations. Sometimes simple exact matching is enough.

#   ==    : returns TRUE for equality of pattern and object, FALSE otherwise.
#           Can only use exact matches, not regexes.
#
# grep()  : returns the index of all matches, can use regexes.
#
# grepl() : returns a logical vector like ==, but can use regexes.
#
# a %in% b: returns a logical vector which is TRUE for an element of a if it
#           exactly matches an element of b.

# Try:
s <- "quid est ergo tempus
      si nemo ex me quaerat scio
      si quaerenti explicare velim nescio"
#     What is therefore time?
#     If no one asks me, I know.
#     If I wish to explain to whoever asks, I do not know
#     (Augustine of Hippo, Confessiones, XI - ca. 397 CE)

v <- unlist(strsplit(s, "\\s+"))
# strsplit(<string>, <regex>) returns a list.
# unlist(<list>) flattens the list into a vector.

#
# The regex "\\s+" explained:
# \s means "whitespace" for the regex engine i.e. " ", "\t", and "\n". However
# inside an R string, \s is not a recognized escape sequence. In order to
# pass the expression "\s" to the regex engine, we nee to pass the literal
# characters "\" and "s". Of course, "s" is just "s", but in order to pass on
# a literal "\" and not have it attach to the next character, we have to
# escape the "\" - and that is done by typing "\\". Therefore "\\s" sends
# the string "\s" to the regex engine, just as "\\t" would specify the string
# "\t" and not a tab character.
# "+" means: one or more matches.
# Therefore "\\s+" means: split the string, consuming all runs of one or more
#                         whitespaces.

v

v == "si"         # using ==
which(v == "si")
v[v == "si"]      # obviously

grep("si", v)     # using grep()
grep("scio", v)
v[grep("scio", v)]    # since "scio" is a regex, not an exact match, partial
                      # matches to elements of v are found. However, this is
                      # still _exact_ matching.
                      #
v[grep("^....$", v)]  # _This_ is ambiguous matching (element contains exactly
                      #                               four of any character.)
                      # N.b - this particular case should better be written:
v[grep("^.{4}$", v)]  # ... or even:
v[nchar(v) == 4]

v %in% c("ut", "ex", "ab", "si")     # using %in%
v[v %in% c("ut", "ex", "ab", "si")]


# ==== Returning matches with regexpr() and regmatches()
# get "o" and the three characters preceding it
pattern <- ".{3}o"
regexpr(pattern, v)      # positions of matches
M <- regexpr(pattern, v) # assign the result to M
regmatches(v, M)         # use regmatches to process M in v

# Note!
# regexpr() captures the FIRST occurrence of the match -
# consider:
patt <- ".a"   # "a" and the preceding character
regmatches(v, regexpr(patt, v)) # ...  but "quaerat" has _two_ "a"!

# You must use gregexpr() to return ALL matches.
patt <- ".a"   # "a" and the preceding character
regmatches(v, gregexpr(patt, v)) # ...  but now the result is a list (obviously)
unlist(regmatches(v, gregexpr(patt, v)))

# Nb. regexepr() and gregexepr() are to each other like sub() and gsub() -
#     the "g-" prefix means "global".


# ==== Capturing matches with regexec()
#
# regexec() captures groups from a string. Such groups are defined by parentheses.
# Example: is there a three-consonant cluster in our string?
# If yes, we want to know the letter preceeding and following it.
patt <- "(.)([bcdfghjklmnpqrstvwxz]{3})(.)"  # Note the parentheses
                                             # that indicate the match
                                             # should be "captured"
grep(patt, s)  # Yes. Ok - what is it?
               # Nb: this is in our source-string s, not the vector v

M <- regexec(patt, s)
regmatches(s, M)                # The result is a list with an object for every
                                # element in s. s has only one ...
regmatches(v, regexec(patt, v)) # ... but v has several objects

regmatches(s, M)[[1]]           # extract from the list
paste(regmatches(s, M)[[1]][2:4], collapse = "-")

# Exercise: is there a three-vowel cluster?


# Unfortunately there is no option to capture multiple matches
# in base R: regexec() lacks a corresponding gregexec()...

# Solution 1 (base R): you can use multiple matches in an sapply()
# statement...

patt <- "(.)(a)(.)"  # "a" environments
sapply(regmatches(s, gregexpr(patt, s))[[1]],
       function(M){regmatches(M, regexec(patt, M))})


# Solution 2 (probably preferred): you can use
# str_match_all() from the very useful library "stringr" ...
if (!require(stringr, quietly=TRUE)) {
    install.packages("stringr")
    library(stringr)
}

str_match_all(s, patt)
str_match_all(v, patt)

# An interesting new alternative/complement to the base R regex libraries is the
# "ore" package that uses the Oniguruma regular expression libraries and
# supports multiple character encodings.
if (!require(ore)) {
    install.packages("ore")
    library(ore)
}

ore.search("o .", s)                 # first match
ore.search("o .", s, all=TRUE)       # all matches
M <- ore.search("o .", s, all=TRUE)
M$nMatches
M$match
str(M)

# According to the author John Clayden, key advantages include:
# - Search results focus around the matched substrings (including
#   parenthesised groups), rather than the locations of matches. This saves
#   extra work with substr(")) or similar to extract the matches themselves.
# - http://rpubs.com/jonclayden/regex-performance Substantially better
#   performance, especially when matching against long strings.
# - Substitutions can be functions as well as strings.
# -  Matches can be efficiently obtained over only part of the strings.
# -  Fewer core functions, with more consistent names.

# === Practice with a sequence
# Get the sequence for MBP1_SACCE:
( s <- refDB$protein[refDB$protein$name == "MBP1_SACCE", "sequence"] )

# ...and paste it into the data box of http://regexpal.com
# Make sure the dialect is set to PCRE, not Javascript.
# Whenever we discuss expressions below, try them out online, and experiment
# with variations.

str_match_all(s, "[DE]{3}")  # clusters of three acidic residues
str_match_all(s, "[HKR]{4}")  # clusters of four basic residues

str_match_all(s, "(([MCILYFW]...){3,})")  # amphipathic helix
ore.search("([MCILYFW]...){3,}", s, all=TRUE)

ore.search("([LIVM].{6}){3,}[LIVM]", s, all=TRUE)  # "leucine" zipper
                                                   # protein-protein interaction
                                                   # motif
# Compare this to MBP1_SPIPU
ore.search("([LIVM].{6}){3,}[LIVM]", myDB$protein$sequence[11], all=TRUE)

# What do you think about the result?
# Doe these results match the pepcoils results?


unlist(strsplit(s, "[KR]"))  # tryptic digest - not.
                             # strsplit() "eats" the regular expression.
                             # We need a look-ahead / look-behind expression
                             # instead, that splits bewteen [KR] and a following
                             # character. The parameter perl = TRUE defines the
                             # proper regex dialect for this:
unlist(strsplit(s, "(?<=[KR])(?=.)", perl = TRUE))

# Refine! Trypsin does not cut before "P"
unlist(strsplit(s, "(?<=[KR])(?=[^P])", perl = TRUE))

# size distribution:
hist(nchar(unlist(strsplit(s, "(?<=[KR])(?=[^P])", perl = TRUE))))

# ==== Unicode?
x <- "ᐅᖃᐅᓯᖅ ᐊᑕᐅᓯᖅ ᓈᒻᒪᔪᐃᑦᑐᖅ"  # Uqauhiq atauhiq naammayuittuq
                                # One language is never enough (Inuktitut)
# Find words ending in ᓯᖅ (iq):

regmatches(x, gregexpr("[^ ]+ᓯᖅ", x))[[1]]
str_match_all(x, "[^ ]+ᓯᖅ")[[1]][, 1]
ore.search("[^ ]+ᓯᖅ", x, all=TRUE)$matches

# Any questions? Curious how to achieve this or that? E-mail the list...




# ==============================================================================
#        PART THREE: SEQUENCE ANALYSIS
# ==============================================================================


if (!require(seqinr, quietly=TRUE)) {
    install.packages("seqinr")
    library(seqinr)
}

help(package = seqinr) # shows the available functions

# Let's try a simple function
?computePI

# This takes as input a vector of upper-case one-letter amino acid codes.
# Let's retrieve two of the MBP1 sequences from our database.

MBP1_SACCE <- myDB$protein[myDB$protein$name == "MBP1_SACCE", "sequence"]
MBP1_SPIPU <- myDB$protein[myDB$protein$name == "MBP1_SPIPU", "sequence"]

# We can use the function strsplit() to split the string
# into single characters

head(unlist(strsplit(MBP1_SACCE, "")), 20) # splitting on the empty spring
                                           # splits into single characters

# Alternatively, seqinr provides
# the function s2c() to convert strings into
# character vectors (and c2s to convert them back).

head(s2c(MBP1_SACCE), 20)

computePI(s2c(MBP1_SACCE))  # isoelectric point
computePI(s2c(MBP1_SPIPU))

pmw(s2c(MBP1_SACCE))        # molecular weight
pmw(s2c(MBP1_SPIPU))        # molecular weight

AAstat(s2c(MBP1_SACCE))     # This also plots the distribution of
AAstat(s2c(MBP1_SPIPU))     # values along the sequence

# ==== Using amino acid data
# A true Labor of Love has gone into the
# compilation of the "aaindex" data:

?aaindex
data(aaindex)  # "attach" the dataset - i.e. make it accessible as an
# R object

length(aaindex)

# Here are all the index descriptions
for (i in 1:length(aaindex)) {
    cat(paste(i, ": ", aaindex[[i]]$D, "\n", sep=""))
}

# ==== A compositional enrichment plot
# Lets use one of the indices to calculate and plot amino-acid
# composition enrichment:
aaindex[[459]]

# Let's construct an enrichment plot to compare one of the amino acid indices
# with the situation in our sequence.

( refData <- aaindex[[459]]$I )      # reference frequencies in %
names(refData) <- a(names(refData))  # change names to single-letter
                                     # code using seqinr's "a()" function
refData


# tabulate our sequence of interest and normalize
( obsData <- table(s2c(MBP1_SACCE)) )        # count occurrences with table()
(obsData <- 100 * (obsData / sum(obsData)) ) # Normalize


logRatio <- numeric() # create an empty vector

# loop over all elements of the reference,
# calculate log-ratios, convert to log2,
# and store them in the vector.
for (aa in names(refData)) {
    logRatio[aa] <- log(obsData[aa] / refData[aa]) / log(2)
}

barplot(logRatio)

# Sort by frequency, descending
logRatio <- sort(logRatio, decreasing = TRUE)

barplot(logRatio)  # If you can't see all of the amino acid letters in the
# x-axis legend, make the plot wider by dragging the
# vertical pane-separator to the left


# label the y-axis
# (see text() for details)
label <- expression(paste(log[2],"( f(obs) / f(ref) )", sep = ""))

barplot(logRatio,
        main = paste("AA composition enrichment"),
        ylab = label,
        cex.names=0.9)



# color the bars by type.
# define colors
chargePlus  <- "#282A8F"
chargeMinus <- "#931523"
hydrophilic <- "#57808F"
hydrophobic <- "#6F7A79"
plain       <- "#FDF7F7"



# Assign the colors to the different amino acid names
barColors <- character()

for (AA in names(logR)) {
    if      (grepl("[HKR]",     AA)) {barColors[AA] <- chargePlus }
    else if (grepl("[DE]",      AA)) {barColors[AA] <- chargeMinus}
    else if (grepl("[NQST]",    AA)) {barColors[AA] <- hydrophilic}
    else if (grepl("[FMILYVW]", AA)) {barColors[AA] <- hydrophobic}
    else                              barColors[AA] <- plain
}

barplot(logRatio,
        main = paste("AA composition enrichment"),
        ylab = label,
        col = barColors,
        cex.names=0.8)


# draw a horizontal line at y = 0
abline(h=0)

# add a legend that indicates what the colours mean
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
                 plain)
)

# == TASK: Interpret this plot. What biological facts does it express about
#          the protein sequence?

# == TASK: I have assembled the plotting code in the function
#          plotAAenrichment(), in the .utilities.R script. Study the code:

plotAAenrichment
#          Use the function to plot the enrichment profile for
#          MBP1_SPIPU. Then, using the back and forward arrows of the Plot Pane,
#          analyze if/how the enrichment profiles are similar or different.
#          What biological facts does this express?

# == TASK: We normally return() values from a function. Here we use invisible()
#          instead. Why? How can we access the returned values?


# ==============================================================================
#        PART FOUR: DOTPLOT AND MDM
# ==============================================================================

# If you come here after restarting RStudio, you need to reload the database.
# load("myDB.01.RData")

# First, we install and load the Biostrings package.
if (!require(Biostrings, quietly=TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
    library(Biostrings)
}


# This is a large collection of tools ...
help(package = "Biostrings")

# ... with mutation matrices and other useful datasets
data(package = "Biostrings")

# Let's load BLOSUM62
data(BLOSUM62)

# ... and see what it contains. (You've seen this before, right?)
BLOSUM62

# We can simply access values via the row/column names to look at the data
# for the questions I asked on the Wiki page:
BLOSUM62["H", "H"]
BLOSUM62["S", "S"]

BLOSUM62["L", "K"]
BLOSUM62["L", "I"]


BLOSUM62["R", "W"]
BLOSUM62["W", "R"]   # the matrix is symmetric!

# Now let's craft code for a dotplot. That's surprisingly simple. We build a
# matrix that has as many rows as one sequence, as many columns as another. Then
# we go through every cell of the matrix and enter the pairscore we encounter
# for the amino acid pair whose position corresponds to the row and column
# index. Finally we visualize the matrix in a plot.
#

# First we fetch our sequences and split them into single characters.
sel <- myDB$protein$name == "MBP1_SACCE"
vSACCE <- s2c(myDB$protein$sequence[sel])

sel <- myDB$protein$name == "MBP1_SPIPU"
vSPIPU <- s2c(myDB$protein$sequence[sel])

# Check that we have two character vectors of the expected length.
str(vSACCE)
str(vSPIPU)

# How do we get the pairscore values? Consider: a single pair of amino acids can
# be obtained from sequence SACCE and YFO eg. from position 13 and 21 ...
vSACCE[13]
vSPIPU[21]

# ... using these as subsetting expressions, we can pull the pairscore
# from the MDM
BLOSUM62[vSACCE[13], vSPIPU[21]]

# First we build an empty matrix that will hold all pairscores ...
dotMat <- matrix(numeric(length(vSACCE) * length(vSPIPU)),
                 nrow = length(vSACCE), ncol = length(vSPIPU))

# ... then we loop over the sequences and store the scores in the matrix.
#
for (i in 1:length(vSACCE)) {
    for (j in 1:length(vSPIPU)) {
        dotMat[i, j] <- BLOSUM62[vSACCE[i], vSPIPU[j]]
    }
}

# Even though this is a large matrix, this does not take much time ...
# Let's have a look at a small block of the values:

dotMat[1:10, 1:10]

# Rows in this matrix correspond to an amino acid from vSACCE, columns in
# the matrix correspond to an amino acid from vSPIPU.

# To plot this, we use the image() function. Here, with default parameters.

image(dotMat)

# Be patient, this takes a few moments to render: more than 800,000 values.
# Nice.
# What do you expect?
# What would similar sequences look like?
# What do you see?

#You migh notice a thin line of yellow along the diagonal, moving approximately
# from bottom left to top right, fading in and out of existence. This is the
# signature of extended sequence similarity.

# Let's magnify this a bit by looking at only the first 200 amino acids ...
image(dotMat[1:200, 1:200])

# ... and, according to our normal writing convention, we would like the
# diagonal to run from top-left to bottom-right since we write from left to
# right and from top to bottom...
image(dotMat[1:200, 1:200], ylim = 1.0:0.0)

# ... and we would like the range of the x- and y- axis to correspond to the
# sequence position ...
image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1))

# ... and labels! Axis labels would be nice ...
image(x = 1:200, y = 1:200,  dotMat[1:200, 1:200], ylim=c(200,1),
      xlab = "MBP1_SPIPU", ylab = "MBP1_SACCE" )

# ... and why don't we have axis-numbers on all four sides? Go, make that right
# too ...
len <- 200
image(x = 1:len, y = 1:len,  dotMat[1:len, 1:len], ylim=c(len,1),
      xlab = "MBP1_SPIPU", ylab = "MBP1_SACCE", axes = FALSE)
box()
axis(1, at = c(1, seq(10, len, by=10)))
axis(2, at = c(1, seq(10, len, by=10)))
axis(3, at = c(1, seq(10, len, by=10)))
axis(4, at = c(1, seq(10, len, by=10)))

# ... you get the idea, we can infinitely customize our plot. However a good way
# to do this is to develop a particular view for, say, a report or publication
# in a script and then put it into a function. I have put a function into the
# utilities file and called it dotPlot2(). Why not dotPlot() ... that's because
# there already is a dotplot function in the seqinr package:

dotPlot(vSACCE, vSPIPU)                                   # seqinr
dotPlot2(vSACCE, vSPIPU, xlab = "SACCE", ylab = "SPIPU")  # Our's

# Which one do you prefer? You can probably see the block patterns that arise
# from segments of repetitive, low complexity sequence. But you probably have to
# look very closely to discern the faint diagonals that correspond to similar
# sequence.


# Let's see if we can enhance the contrast between distributed noise and the
# actual alignment of conserved residues. We can filter the dot matrix with a
# pattern that enhances diagonally repeated values. Every value in the matrix
# will be replaced by a weighted average of its neighborhood. Here is  a
# diagonal-filter:

myFilter <- matrix(numeric(25), nrow = 5)
myFilter[1, ] <- c( 1, 0, 0, 0, 0)
myFilter[2, ] <- c( 0, 1, 0, 0, 0)
myFilter[3, ] <- c( 0, 0, 1, 0, 0)
myFilter[4, ] <- c( 0, 0, 0, 1, 0)
myFilter[5, ] <- c( 0, 0, 0, 0, 1)

# I have added the option to read such filters (or others that you could define on your own) as a parameter of the function.

dotPlot2(vSACCE, vSPIPU, xlab = "SACCE", ylab = "SPIPU", f = myFilter)

# I think the result shows quite nicely how the two sequences are globally
# related and where the regions of sequence similarity are. Play with this a bit
# ...  Can you come up with a better filter? If so, eMail us.



# ==============================================================================
#        PART FIVE: BIOSTRINGS PAIRWISE ALIGNMENT
# ==============================================================================

# Biostrings is one of the basic packages that the Bioconductor software
# landscape builds on. It stores sequences in "AAstring" objects and these are
# complex software structures that are designed to be able to handle
# genome-scale sequences. Biostrings functions - such as the alignment functions
# - expect their input to be Biostrings objects.
?AAString

AAString("ACDE")
s <- AAString("ACDE")
str(s)
# See: it's complicated. This is an "S4" object. Bioconductor uses these objects
# almost exclusively, but we will not be poking around in their internals. Just
# this: how do we get the sequence back out of an AAString object? The help page
# for XString - the parent "class" of AAStrings - mentions the  alternatives:

as.character(s)  # the base R version
toString(s)      # using the Biostrings function toString()

# While we need to remember to convert our sequences from the character vectors
# that we store in our database, to AAStrings that we can align, the alignment
# itself is really straightforward. The pairwiseAlignment() function was written
# to behave exactly like the functions you encountered on the EMBOSS server.

# First: make AAString objects ...
sel <- myDB$protein$name == "MBP1_SACCE"
aaMBP1_SACCE <- AAString(myDB$protein$sequence[sel])

sel <- myDB$protein$name == "MBP1_SPIPU"
aaMBP1_SPIPU <- AAString(myDB$protein$sequence[sel])


?pairwiseAlignment

# ... and align.
# Global optimal alignment with end-gap penalties is the default.
# (like EMBOSS needle)
ali1 <-  pairwiseAlignment(
    aaMBP1_SACCE,
    aaMBP1_SPIPU,
    substitutionMatrix = "BLOSUM62",
    gapOpening = 10,
    gapExtension = 0.5)

str(ali1)  # Did you think the AAString object was complicated ?

# This is a Biostrings alignment object. But we can use Biostrings functions to
# tame it:
ali1
writePairwiseAlignments(ali1)   # That should look familiar

# And we can make the internal structure work for us  (@ is for classes as
# $ is for lists ...)
str(ali1@pattern)
ali1@pattern
ali1@pattern@range
ali1@pattern@indel
ali1@pattern@mismatch

# or work with "normal" R functions
# the alignment length
nchar(ali1@pattern)

# the number of identities
sum(s2c(as.character(ali1@pattern)) ==
        s2c(as.character(ali1@subject)))

# ... e.g. to calculate the percentage of identities
100 *
    sum(s2c(as.character(ali1@pattern)) ==
            s2c(as.character(ali1@subject))) /
    nchar(ali1@pattern)
# ... which should be the same as reported in the writePairwiseAlignments()
# output. Awkward to type? Then it calls for a function:
#
percentID <- function(al) {
    # returns the percent-identity of a Biostrings alignment object
    return(100 *
               sum(s2c(as.character(al@pattern)) ==
                       s2c(as.character(al@subject))) /
               nchar(al@pattern))
}

percentID(ali1)

# Compare with local optimal alignment (like EMBOSS Water)
ali2 <-  pairwiseAlignment(
    aaMBP1_SACCE,
    aaMBP1_SPIPU,
    type = "local",
    substitutionMatrix = "BLOSUM62",
    gapOpening = 50,
    gapExtension = 10)

writePairwiseAlignments(ali2)   # This has probably only aligned the N-terminal
# DNA binding domain - but that one has quite
# high sequence identity:
percentID(ali2)

# == TASK: ==

# Compare the two alignments. I have weighted the local alignment heavily
# towards an ungapped alignment by setting very high gap penalties. Try changing
# the gap penalties and see what happens: how does the number of indels change,
# how does the length of indels change...

# Fine. Please return to the Wiki to study BLAST alignment...


# ==============================================================================
#        PART SIX: APSES DOMAIN ANNOTATION BY ALIGNMENT
# ==============================================================================

# In this section we define the SPIPU APSES domain sequence by performing a global,
# optimal sequence alignment of the yeast domain with the full length protein
# sequence of MBP1_SPIPU.
#

# I have annotated the yeast APSES domain as a proteinAnnotation in the
# database. To view the annotation, we can retrieve it via the proteinID and
# featureID. Here is the yeast protein ID:
myDB$protein$ID[myDB$protein$name == "MBP1_SACCE"]

# ... assign it for convenience:
proID <- myDB$protein$ID[myDB$protein$name == "MBP1_SACCE"]

# ... and if you look at the feature table, you can identify the feature ID
myDB$feature[ , c("ID", "name", "description")]
myDB$feature$ID[myDB$feature$name == "APSES fold"]

# ... assign it for convenience:
ftrID <- myDB$feature$ID[myDB$feature$name == "APSES fold"]

# ... and with the two annotations we can pull the entry from the protein
# annotation table
myDB$proteinAnnotation[myDB$proteinAnnotation$protein.ID == proID &
                       myDB$proteinAnnotation$feature.ID == ftrID, ]

myDB$proteinAnnotation$ID[myDB$proteinAnnotation$protein.ID == proID &
                          myDB$proteinAnnotation$feature.ID == ftrID]

# ... assign it for convenience:
fanID <- myDB$proteinAnnotation$ID[myDB$proteinAnnotation$protein.ID == proID &
                                   myDB$proteinAnnotation$feature.ID == ftrID]

# The annotation record contains the start and end coordinates which we can use
# to define the APSES domain sequence with a substr() expression.
substr(myDB$protein$sequence[myDB$protein$ID == proID],
       myDB$proteinAnnotation$start[myDB$proteinAnnotation$ID == fanID],
       myDB$proteinAnnotation$end[myDB$proteinAnnotation$ID == fanID])

# Lots of code. But don't get lost. Let's recapitulate what we have done: we
# have selected from the sequence column of the protein table the sequence whose
# name is "MBP1_SACCE", and selected from the proteinAnnotation table the start
# and end coordinates of the annotation that joins an "APSES fold" feature with
# the sequence, and used the start and end coordinates to extract a substring.
# The expressions get lengthy, but it's not hard to wrap all of this into a
# function so that we only need to define name and feature.

dbGetFeatureSequence
dbGetFeatureSequence(myDB, "MBP1_SACCE", "APSES fold")


# Let's convert this to an AAstring and assign it:
aaMB1_SACCE_APSES <- AAString(dbGetFeatureSequence(myDB,
                                                   "MBP1_SACCE",
                                                   "APSES fold"))

# To align, we need the YFO sequence. Here is it's definition again, just
# in case ...

sel <- myDB$protein$name == "MBP1_SPIPU"
aaMBP1_SPIPU <- AAString(myDB$protein$sequence[sel])

# Now let's align these two sequences of very different length without end-gap
# penalties using the "overlap" type. "overlap" turns the
# end-gap penalties off and that is crucially important since
# the sequences have very different length.

aliApses <-  pairwiseAlignment(
    aaMB1_SACCE_APSES,
    aaMBP1_SPIPU,
    type = "overlap",
    substitutionMatrix = "BLOSUM62",
    gapOpening = 10,
    gapExtension = 0.5)

# Inspect the result. The aligned sequences should be clearly
# homologous, and have (almost) no indels. The entire "pattern"
# sequence from QIYSAR ... to ... KPLFDF  should be matched
# with the "query". Is this correct?
writePairwiseAlignments(aliApses)

# If this is correct, you can extract the matched sequence from
# the alignment object. The syntax is a bit different from what
# you have seen before: this is an "S4 object", not a list. No
# worries: as.character() returns a normal string.
as.character(aliApses@subject)

# Now, what are the aligned start and end coordinates? You can read them from
# the output of writePairwiseAlignments(), or you can get them from the range of
# the match.

str(aliApses@subject@range)

# start is:
aliApses@subject@range@start

# ... and end is:
aliApses@subject@range@start + aliApses@subject@range@width - 1

# Since we have this section defined now, we can create a feature annotation
# right away and store it in myDB.
myProteinID <- "my_pro_1"
myFeatureID <- "ref_ftr_1"
myStart <- 22
myEnd   <- 120

# ==== create the proteinAnnotation entry
panRow <- data.frame(ID = dbCreateID(myDB$proteinAnnotation$ID),
                     protein.ID = myProteinID,
                     feature.ID = myFeatureID,
                     start = myStart,
                     end = myEnd,
                     stringsAsFactors = FALSE)
myDB$proteinAnnotation <- rbind(myDB$proteinAnnotation, panRow)

# == check that this was successful and has the right data
myDB$proteinAnnotation[nrow(myDB$proteinAnnotation), ]

# ===== end code template ===========================================

# ... continue here.
# I expect that a correct result would look something like
#           ID protein.ID feature.ID start end
# 105 my_fan_1   my_pro_1  ref_ftr_1    22 120

# If you made a mistake, simply overwrite the current version of myDB by loading
# your saved, good version:  load("myDB.01.RData") and correct your mistake.

# If this is correct, save it
save(myDB, file = "myDB.02.RData")  # Note that it gets a new version number!

# Done with this part. Copy the sequence of the APSES domain of MBP1_SPIPU
as.character(aliApses@subject)
# ... you need it for the reverse BLAST search, and return to the Wiki page.


# ==============================================================================
#        PART FIVE: ADD PSI BLAST RESULTS
# ==============================================================================

# Once we have discovered APSES domain proteins in our species, we should enter them
# into the database. Execute the code below to add the additional SACCE proteins ...

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SWI4_SACCE",
              RefSeqID = "NP_011036",
              UniProtID = "P25302",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
1 mpfdvlisnq kdntnhqnit pisksvllap hsnhpvieia tysetdvyec yirgfetkiv
61 mrrtkddwin itqvfkiaqf sktkrtkile kesndmqhek vqggygrfqg twipldsakf
121 lvnkyeiidp vvnsiltfqf dpnnpppkrs knsilrktsp gtkitspssy nktprkknss
181 sstsatttaa nkkgkknasi nqpnpsplqn lvfqtpqqfq vnssmnimnn ndnhttmnfn
241 ndtrhnlinn isnnsnqsti iqqqksihen sfnnnysatq kplqffpipt nlqnknvaln
301 npnnndsnsy shnidnvins snnnnngnnn nliivpdgpm qsqqqqqhhh eyltnnfnhs
361 mmdsitngns kkrrkklnqs neqqfynqqe kiqrhfklmk qpllwqsfqn pndhhneycd
421 sngsnnnnnt vasngssiev fssnendnsm nmssrsmtpf sagntssqnk lenkmtdqey
481 kqtiltilss erssdvdqal latlypapkn fninfeiddq ghtplhwata maniplikml
541 itlnanalqc nklgfncitk sifynncyke nafdeiisil kiclitpdvn grlpfhylie
601 lsvnksknpm iiksymdsii lslgqqdynl lkiclnyqdn igntplhlsa lnlnfevynr
661 lvylgastdi lnldnespas imnkfntpag gsnsrnnntk adrklarnlp qknyyqqqqq
721 qqqpqnnvki pkiiktqhpd kedstadvni aktdsevnes qylhsnqpns tnmntimedl
781 sninsfvtss vikdikstps kilenspily rrrsqsisde kekakdnenq vekkkdplns
841 vktampsles pssllpiqms plgkyskpls qqinklntkv sslqrimgee iknldnevve
901 tessisnnkk rlitiahqie dafdsvsnkt pinsisdlqs riketsskln sekqnfiqsl
961 eksqalklat ivqdeeskvd mntnssshpe kqedeepipk stsetsspkn tkadakfsnt
1021 vqesydvnet lrlateltil qfkrrmttlk iseakskins svkldkyrnl igitienids
1081 klddiekdlr ana"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "PHD1_SACCE",
              RefSeqID = "NP_012881",
              UniProtID = "P36093",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
1 myhvpemrlh yplvntqsna aitptrsydn tlpsfnelsh qstinlpfvq retpnayanv
61 aqlatsptqa ksgyycryya vpfptypqqp qspyqqavlp yatipnsnfq pssfpvmavm
121 ppevqfdgsf lntlhphtel ppiiqntndt svarpnnlks iaaasptvta ttrtpgvsst
181 svlkprvitt mwedenticy qveangisvv rradnnming tkllnvtkmt rgrrdgilrs
241 ekvrevvkig smhlkgvwip ferayilaqr eqildhlypl fvkdiesivd arkpsnkasl
301 tpksspapik qepsdnkhei ateikpksid alsngastqg agelphlkin hidteaqtsr
361 aknels"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SOK2_SACCE",
              RefSeqID = "NP_013729",
              UniProtID = "P53438",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
1 mpignpintn diksnrmrqe snmsavsnse stigqstqqq qqqqqylgqs vqplmpvsyq
61 yvvpeqwpyp qyyqqpqsqs qqqlqsqpqm yqvqesfqss gsdsnasnpp stsvgvpsna
121 tatalpngsa ittkksnnst nisnnvpyyy yfpqmqaqqs maysypqayy yypangdgtt
181 ngatpsvtsn qvqnpnlekt ystfeqqqqh qqqqqlqaqt ypaqppkign afskfsksgp
241 psdsssgsms pnsnrtsrns nsisslaqqp pmsnypqpst yqypgfhkts sipnshspip
301 prslttptqg ptsqngplsy nlpqvgllpp qqqqqvsply dgnsitppvk pstdqetylt
361 anrhgvsdqq ydsmaktmns fqtttirhpm pliattnatg sntsgtsasi irprvtttmw
421 edektlcyqv eangisvvrr adndmvngtk llnvtkmtrg rrdgilkaek irhvvkigsm
481 hlkgvwipfe ralaiaqrek iadylyplfi rdiqsvlkqn npsndsssss sstgiksisp
541 rtyyqpinny qnpngpsnis aaqltyssmn lnnkiipnns ipavstiaag ekplkkctmp
601 nsnqleghti tnlqtlsatm pmkqqlmgni asplsyprna tmnsastlgi tpadskpltp
661 sptttntnqs sesnvgsiht gitlprvese sashskwske adsgntvpdn qtlkeprssq
721 lpisaltstd tdkiktstsd eatqpnepse aepvkesess ksqvdgagdv sneeiaaddt
781 kkqek"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "XBP1_SACCE",
              RefSeqID = "NP_012165",
              UniProtID = "P40489",
              taxonomy.ID = as.integer(4932),
              sequence = dbSanitizeSequence("
1 mkypafsins dtvhltdnpl ddyqrlylvs vldrdsppas fsaglnirkv nykssiaaqf
61 thpnfiisar dagngeeaaa qnvlncfeyq fpnlqtiqsl vheqtllsql assatphsal
121 hlhdknilmg kiilpsrsnk tpvsasptkq ekkalstasr enatssltkn qqfkltkmdh
181 nlindklinp nncviwshds gyvfmtgiwr lyqdvmkgli nlprgdsvst sqqqffckae
241 fekilsfcfy nhssftsees ssvllsssts sppkrrtstg stfldanass sstsstqann
301 yidfhwnnik pelrdlicqs ykdflinelg pdqidlpnln panftkrirg gyikiqgtwl
361 pmeisrllcl rfcfpiryfl vpifgpdfpk dceswylahq nvtfassttg agaataataa
421 antstnftst avarprqkpr prprqrstsm shskaqklvi edalpsfdsf venlglssnd
481 knfikknskr qksstytsqt sspigprdpt vqilsnlasf ynthghrysy pgniyipqqr
541 yslpppnqls spqrqlnyty dhihpvpsqy qsprhynvps spiapapptf pqpygddhyh
601 flkyasevyk qqnqrpahnt ntnmdtsfsp rannslnnfk fktnskq"),
              stringsAsFactors = FALSE))

# ... then execute the code below to add the additional SPIPU proteins.

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SPPG_08986",
              RefSeqID = "XP_016611096",
              UniProtID = "A0A0L0HPN9",
              taxonomy.ID = as.integer(645134),
              sequence = dbSanitizeSequence("
        1 mirahdmafe ggaqypslsh ttvgssgaga iatpangsls dtlagglstp lpplssshdf
       61 lqsyfgnssl nqmsgssgpp pflplqipal srltegnsil gvsprniynq ppspkrlpee
      121 eseedeghad ggssttstpa palngegeiy aavysgipvy emmcnsvavm rrksdgflna
      181 tqilkvagvd kgkrtkildk evmigehekv qggygkyqgt wipfergvql aamygvdgfl
      241 rpllefelpa pgradqtptk eqvmaanrdl lkrsnsssfn ksksrshgdd gtdtaqrrrr
      301 agapstrask rlmeiessdi hddessdnpl aspgplgnya kkrprtdaya esvietnaek
      361 yramlmamfv hedplyipdm lsgptlpndl dldivideqs htslhwaaal arinvvrvll
      421 qkgadirsvn ndgetalira vrvtnnydnq tfpelldllh stlslvdakn rtvlhhiaat
      481 aglegrvaas ryylecllew varhggnfss lvdiqdsagd talnvaarig nrnlveqlid
      541 vgadgtienr aglrpsdfgf edvlgssrgd ggvegtgdtg eaekiegtsr vvfpsirtee
      601 edaalalsav asatkgreia savqqmvdem sctfsaemkv kadqltearn qlrqvtkela
      661 dvrkqnhilr kenqmlpemv qqiknlercl geemakcnta tihdqrstds tpggmvtggi
      721 tetkqtpqdp ngelaalrrt lqekealers lrseiirlks tsgkselgck kiiaaccnvp
      781 vesvdellnp llqavesdke vdmsavagfm mdvkrreggr c"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SPPG_02210",
              RefSeqID = "XP_016611186",
              UniProtID = "A0A0L0HPY0",
              taxonomy.ID = as.integer(645134),
              sequence = dbSanitizeSequence("
        1 mstslyepnf laslkddaed avdpaelvlp tpppppppsy ptalvppvgp pyqpwedpry
       61 algpmphykp vyvqqhqahv pqpappapey ytppaaspss tastlqnilp qaptssymmp
      121 pqhlmhqqql hyysmqaqpp vnvpvfrpsn fpipalrrnn sshmgtpakr ptpslpprmp
      181 vparkglyea tysnvpvyem aspngvgvmr rrhdswmnat hilkaagmek srrtkvlere
      241 ihqgehekiq ggygkyqgtw iplyrarela aeygldealr dildlpe"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SPPG_09382",
              RefSeqID = "XP_016606105",
              UniProtID = "A0A0L0H980",
              taxonomy.ID = as.integer(645134),
              sequence = dbSanitizeSequence("
        1 mapairaasy sgvrvyeiqc rgvpvmrrva dayinatqil raaglpkpqr tkilerdvtg
       61 gvhekvqggy agfqgtwipl esartlaaah gveqdlvelf dykpqpgdae ppakekrkrp
      121 appkeepkek wssssfvsss dgegtrpfst fvkqvdrgfn stfewsggde qqregsledd
      181 ietpearrtt rmalrppgtr lkrrpvypge deeegaylsp kkrrtespfs haaspmgpst
      241 vdvtksvaaf tqpvrapsfv rdktpphagk rcegcgvtst pqwrrgpsgk rtlcnacgvk
      301 wsfgrlntkg ksqddynsdd sdlmdnemdl didgpdlvpe sspmrpvpla ekttrdlemq
      361 vkhlkrqlre sersrkrlrk ileeavtdds emdrnfrqvv nkctiraavy tppeyaptai
      421 dymsgdddsf sdedehealv vtrflravkq nrarmeynis rmhlnkdwav napnearapd
      481 aimv"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SPPG_02885",
              RefSeqID = "XP_016610458",
              UniProtID = "A0A0L0HNM3",
              taxonomy.ID = as.integer(645134),
              sequence = dbSanitizeSequence("
        1 mlstqlpdrn llengpeded sstilaskpe dmldephpkp rdpasvndyn esnrndstlg
       61 lqlpssldll ahaatesftr krkmsesdtv dkekdgsvdd skrcradqvs elhsqdptll
      121 fnlltqersn naslikdlel lhtshaaknl hitqleavne elvnllverd rriedlvral
      181 vggagsssaa aprsephitt sitlkpskti ianvpttptd slpptpdtti ispqqppals
      241 tanwggipvy qlhlasttim rrktdswinv thvlkaanip rgsrlrilek emhtgeheki
      301 qggygkyqgt wvplhraykl aerygllqtl spifnvecfi re"),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SPPG_00954",
              RefSeqID = "XP_016611510",
              UniProtID = "A0A0L0HPW3",
              taxonomy.ID = as.integer(645134),
              sequence = dbSanitizeSequence("
        1 mdgdiaqqvq kfadlcneaf ggqtiegrln fvrpedtrmh fqlpasfqls tdkmimkrlk
       61 vaqsrnkpye rlpvkssvdg gsqnrkrrms tgsaslgtsr yndyqafqic vngtpvmrrk
      121 sdgwinlthi ckaagmpkgq rsrilskelr hghehekvqg gaasfqgtwv padhakyfai
      181 kngiyeqiqp lfahqdaeft nvmppqhpta vekqlpmltp adslprtpgn pkalnfvels
      241 hvsqmpamny qamnyqcgfh lpypqqqhei gtghstsvmp tpdaddvtkf llsasahplv
      301 inam                                            "),
              stringsAsFactors = FALSE))

myDB$protein <-
    rbind(myDB$protein,
          data.frame(
              ID = dbCreateID(myDB$protein$ID),
              name = "SPPG_07480",
              RefSeqID = "XP_016605125",
              UniProtID = "A0A0L0H6I9",
              taxonomy.ID = as.integer(645134),
              sequence = dbSanitizeSequence("
        1 mvepsscafg ksdtdsnimp sltnslssaf vlkpfpaask tvltyqhfrv emhsddakdk
       61 vflaivkall rqgnrpctpk elsnlilrhk lttlggatpy atvssrisqh fkratelsrp
      121 pllgrrslde rnsrrlvyyv dqvgvpvkpl delahsesde smdvgdhlsd sfvsslassp
      181 vlapgsidil sgahgdtgis sslkegaggg adgpatplrv strvkrrrsp ysppddqqsk
      241 rsrsaslasa vtpivspssg tgtddmdkpv ilvedvelvd lstpaasrrs svvdasvpev
      301 ganiepdsik nedpttdseq dfhedmmkvd geygsgrsyq parlipstps tyrqrvssde
      361 fllkpitspq ipsfqsslps klspsrphlg fepadhtips pffdpevdpl gveahesrag
      421 pkaenysdlr mfdvgiheel dyhefhhped msvseldlll redegswgvl ggvtpgqspr
      481 krpaeslapt ptrmrdsvvt qrqdagsesv qvakkrkvqp acegsvndtg lkndksgegd
      541 gdasgqstkr ppvtaeshas tsvsqqsast ipenkqdssp rqnsteqsrs dnvedaatdh
      601 lstticnqav shtpspprtp rndtpsdhsl psrshsalsv asyrcngdtl elyekahdpi
      661 strvskiksa prmelvlrtl hvdslgrail kppsdnepdr kvgirmrgyv rarelfevap
      721 sqdvsvptsp evkvededge didmedesvw vayteevrrr avgprtgdes avvrivgena
      781 pgsetghvlv vlkgsdvegd glvpvecggv wiptsearrl atklgvsaqf nelldmgdes
      841 nngeievlps fstltspklp plpstptrps arslpspast pipsalsskv pdaslgpaml
      901 eddmevdegh dflvfegdee wshdessavt gsegpesnvv sngameklpv tidhstleli
      961 apvqpppnst lplymatidg vlvyitwltk gsdgsiaptp areisdples spagseesga
      1021 atppsapaap tsgptstppp vppppgsssa vpllrrvdnn mvnatlllha gglvtdkers
      1081 ivlslergra rcrkkgsgly gtwiplarar hlarsfclen rlslflsdgf grvvfgvdgv
      1141 ngargkakla nddasgsegt pvaadgsgss vdagseagds sfvdgaasnl tpaaaalaas
      1201 kaamlgltsl ppnvvsglas igagggaagi rglmanrgra rgrpplkhhp isgrgfmnvg
      1261 rlgatlgssp sytpsgtgtn vsttptssam aatqaaiaal sklggngplt aatlasalaa
      1321 lnqvnqqnit agkvttpgss pstptkpaya tgtspspatt ttatssstap aasgiptpas
      1381 ilaaftnmlk ngfpggwkpa gapsintstt sttsappssa sssitfpssl kvdpkpgvss
      1441 splsttpvtt qgvvaralaa agagrgggls gtglsnpvsi qtnalaglpn lgnlggtsfs
      1501 llpftmqnsg nnlslatstt stsatsstgi siggidgptp alplltpptv piiaiddidg
      1561 gdanmpiasp aalaaaaaaa gvtlpvgvey veeeeeeede eeveevteeq geaedenege
      1621 taqggggksr pgakangttk ieqenesdes depdnqdnda ddddddsdie ivkssggksa
      1681 psrsggkirp raavpptsrr apstsgksrg gstgmatrst rsastpntka prrpppapkk
      1741 grvksksrrg attiinhsnk hstvvnngdv gseegsdfei dvvdgdgvdd fr"),
              stringsAsFactors = FALSE))

# TASK: ==
# Confirm that the proteins have been added to myDB with all of their
# required information. What's the easiest way to do that?

# To recapitulate: in addition to the ten original proteins in refDB, myDB
# should now contain MBP1_SPIPU, four more SACCE proteins and six more SPIPU
# proteins.

# If all is correct, save the new version of myDB as version 03: "myDB.03.RData"
save(myDB, file = "myDB.03.RData")



# ==============================================================================
#        PART SIX: SMART DOMAIN ANNOTATIONS
# ==============================================================================

# Plot domain annotations as colored rectangles on a sequence.

# Step one: print the sequence
chunkSeq(myDB$protein$sequence[myDB$protein$name == "MBP1_SPIPU"])

# ... copy, and paste it into the SMART form.

# Step two: enter your domain annotations as features into the database.
#
# load("myDB.03.RData") in case you have restarted the RStudio session.

# == Update myDB

# Every annotated feature requires its own entry in the database. You have added
# the feature for the "APSES fold" before, but copying that code for each feature is clunky and violates our DRY (Do not Repeat Yourself) principle. I
# have constructed the appropriate dbAddFeature function.

# Here is how you call it to add the KilA-N domain annotation for MBP1_SPIPU. You need to know the protein ID:
myDB$protein$ID[myDB$protein$name == "MBP1_SPIPU"]
#... and the feture IDs:
myDB$feature[ , c("ID", "name", "description")]

# ... and this is the function call.
myDB$proteinAnnotation <- dbAnnotateFeature(db = myDB,
                                            pID = "my_pro_1",
                                            fID = "ref_ftr_2",
                                            Start = 40,
                                            End = 123)

# Add every SMART annotated feaure of MBP1_SPIPU, for which there is an entry in
# the feature table of the data model, to the database. Don't do this on the fly
# in the console, but edit your code in myScript.R If you make mistakes, just
# reload the latest version (probably "myDB.03.RData"), then run your corrected
# annotation script again. I see SMART returning:

# In confidently predicted domains:
#  -  a KilAN domain (already added above)
#  -  eight low-complexity regions
#  -  two Ankyrin domains
#  -  two coiled coil domains
# In homologues of know structure:
#  -  one Swi6 fold (annoted as PDB:1SW6|B)
#  In "features NOT shown (overlap only):
#  -  ...

# When you are done adding the additional thirteen features, execute ...
myDB$proteinAnnotation
# ... to confirm that everything is correct.
#
# Once you are sure your annotations are correct, save the database again.
save(myDB, file = "myDB.04.RData") # save the new version


# Then let's plot and compare annotations.
#
# We need a small utility function that draws the annotation boxes on a
# representation of sequence. It will accept the left and right boundaries, the
# height and the color of the box and plot it using R's rect() function.

drawBox <- function(xLeft, xRight, y, colour) {
    # Draw a box from xLeft to xRight at y, filled with colour
    rect(xLeft, (y - 0.1), xRight, (y + 0.1),
         border = "black", col = colour)
}

# test this:
plot(c(-1.5, 1.5), c(0, 0), type = "l")
drawBox(-1, 1, 0.0, "peachpuff")

# Next, we define a function to plot annotations for one protein: the name of
# the protein, a horizontal grey line for its length, and all of its features.

plotProtein <- function(DB, ID, y) {
    # DB: protein database, probably you want myDB
    # ID: the ID of the protein to plot.
    # y: where to draw the plot
    #
    # Define colors: we create a vector of color values, one for
    # each feature, and we give it names of the feature ID. Then we
    # can easily get the color value from the feature name.
    # A: make a vector of color values. The syntax may appear unusual -
    #    colorRampPalette() returns a function, and we simply append
    #    the parameter (number-of-features) without assigning the function
    #    to its own variable name.
    ftrCol <- colorRampPalette(c("#f2003c", "#F0A200", "#f0ea00",
                                 "#62C923", "#0A9A9B", "#1958C3",
                                 "#8000D3", "#D0007F"),
                               space="Lab",
                               interpolate="linear")(nrow(DB$feature))
    # B: Features may overlap, so we make the colors transparent by setting
    #    their "alpha channel" to 1/2  (hex: 7F)
    ftrCol <- paste(ftrCol, "7F", sep = "")
    # C: we asssign names
    names(ftrCol) <- DB$feature$ID
    # E.g. color for the third feature: ftrCol[ DB$feature$ID[3] ]

    # find the row-index of the protein ID in the protein table of DB
    iProtein <- which(DB$protein$ID == ID)

    # write the name of the protein
    text(-30, y, adj=1, labels=DB$protein$name[iProtein], cex=0.75 )

    #draw a line from 0 to nchar(sequence-of-the-protein)
    lines(c(0, nchar(DB$protein$sequence[iProtein])), c(y, y),
          lwd=3, col="#999999")

    # get the rows of feature annotations for the protein
    iFtr <- which(DB$proteinAnnotation$protein.ID == ID)

    # draw a colored box for each feature
    for (i in iFtr) {
        drawBox(DB$proteinAnnotation$start[i],
                DB$proteinAnnotation$end[i],
                y,
                ftrCol[ DB$proteinAnnotation$feature.ID[i] ])
    }
}

# Plot each annotated protein:
# Get the rows of all unique annotated protein IDs in the protein table
iRows <- which(myDB$protein$ID %in% unique(myDB$proteinAnnotation$protein.ID))
# define the size of the plot-frame to accomodate all proteins
yMax <- length(iRows) * 1.1
xMax <- max(nchar(myDB$protein$sequence[iRows])) * 1.1  # longest sequence

# plot an empty frame
plot(1,1, xlim=c(-200, xMax), ylim=c(0, yMax),
     type="n", axes=FALSE, bty="n", xlab="sequence position", ylab="")
axis(1, at = seq(0, xMax, by = 100))

# Finally, iterate over all proteins and call plotProtein()
for (i in 1:length(iRows)) {
    plotProtein(myDB, myDB$protein$ID[iRows[i]], i)
}

# The plot shows clearly what is variable and what is constant about the
# annotations in a group of related proteins. I may ask you to print and submit
# the plot for credit later in class.
#

# ==============================================================================
#        PART SEVEN: MULTIPLE SEQUENCE ALIGNMENT
# ==============================================================================

# We will compute a multiple sequence alignment using the "muscle" algorithm
# which is available throught the Bioconductor msa package.

if (!require(Biostrings, quietly=TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
}
data(BLOSUM62)

if (!require(msa, quietly=TRUE)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("msa")
    library(msa)
}

# If the installation asks you if you want to update older packages, always
# answer "a" for "all" unless you have an important reason not to. But if the
# installer asks you whether you want to compile from source, answer"n" for "no"
# unless you need new functionality in a particular bleeding-edge version of a
# package.

help(package = "msa")

# We have used biostrings' AAString() function before; for multiple
# alignments we need an AAStringSet(). We can simply feed it
# a vector of sequences:

# Let's make a shorthand for selection of Mbp1 proteins from our database: a
# vector of indices for all of the rows in the protein table that contain
# proteins whose name begins with MBP1.
( iMbp1Proteins <- grep("^MBP1_", myDB$protein$name) )

# Next we create an AAStringSet for all of those proteins
( seqSet <- AAStringSet(myDB$protein$sequence[iMbp1Proteins]) )

# ... and align (which is very simple) ...
msaMuscle(
    seqSet,
    order = "aligned")

# ... but to help us make sense of the alignment we need to add the names for
# the sequences. Names for a seqSet object are held in the ranges slot...

seqSet@ranges@NAMES <- myDB$protein$name[iMbp1Proteins]

seqSet

# This little step of adding names is actually really very important. The
# aligned sequences are meaningless strings of characters unless we can easily
# identify their biological relationships. Creating MSAs that are only
# identified by e.g. their RefSeq ids is a type of cargo-cult bioinformatics
# that we encounter a lot. The point of the alignment is not to create it, but
# to interpret it!


# Let's align again, and assign the result ...
( msa1 <-  msaMuscle(
    seqSet,
    order = "aligned") )


# ... or, to see the whole thing (cf. ?MsaAAMultipleAlignment ... print method):
print(msa1, show=c("alignment", "complete"), showConsensus=FALSE)


# You see that the alignment object has sequence strings with hyphens as
# indel-characters. The names are printed to the console. And you also see that
# the order has not been preserved, but the most similar sequences are now
# adjacent to each other.

# Lets write the alignment to one of the common file formats: a multi-fasta
# file.

# Why oh why does the msa package not have a function to do this !!!
# Like, seriously ...

# ==== writeMFA

# Here is our own function to write an AAStringSet object to a multi-FASTA file
# or to block-wise formatted output (source in the .utilities.R script).

# Study:
writeSeqSet

# confirm that the function works
writeSeqSet(seqSet)
writeSeqSet(msa1, format = "ali")

# We can use this function to write the raw and aligned sequences to file like
# so:
writeSeqSet(seqSet, file = "APSES_proteins.mfa")

# writeSeqSet(msa1, file = "APSES_proteins_muscle.mfa")

# ... do that, and continue with the Wiki page section on computing alignments
# online at the EBI.

# == Task:
# Print the output of print(msa1) on a sheet of paper, I may ask you to hand it
# in at a later time for credit.
#



# ==============================================================================
#      PART EIGHT: DOWNLOADING DATA FROM THE WEB
# ==============================================================================

# In previous work, we added sequences to our database by hand, but since the
# information should be cross-referenced and available based on a protein's
# RefSeq ID, we should really have a function that automates this process. It is
# far too easy to make mistakes and enter erroneous information otherwise.

# To begin, we load some libraries with functions
# we need...

# httr sends and receives information via the http protocol, just like a Web
# browser.
if (!require(httr, quietly=TRUE)) {
    install.packages("httr")
    library(httr)
}

# NCBI's eUtils send information in XML format; we need to be able to parse XML.
if (!require(XML, quietly=TRUE)) {
    install.packages("XML")
    library(XML)
}

# stringr has a number of useful utility functions to work with strings. E.g. a
# function that strips leading and trailing whitespace from strings.
if (!require(stringr, quietly=TRUE)) {
    install.packages("stringr")
    library(stringr)
}


# We will walk through the process of fetching data from the Web interfaces of
# UniProt and NCBI with the refSeqID of yeast Mbp1
refSeqID <- "NP_010227"


# === UniProt ID mapping service

# The UniProt ID mapping service supports a "RESTful API": responses can be
# obtained simply via a Web-browser's request. Such requests are commonly sent
# via the GET or POST verbs that a Webserver responds to, when a client asks for
# data. GET requests are visible in the URL of the request; POST requests are
# not directly visible, they are commonly used to send the contents of forms, or
# when transmitting larger, complex data items. The UniProt ID mapping sevice
# can accept long lists of IDs, thus using the POST mechanism makes sense.

# R has a POST() function as part of the httr package.It's very straightforward
# to use: just define the URL of the server and send a list of items as the body
# of the request.

URL <- "http://www.uniprot.org/mapping/" # uniProt ID mapping service
response <- POST(URL,
                 body = list(from = "P_REFSEQ_AC",
                             to = "ACC",
                             format = "tab",
                             query = refSeqID))

response

# If the query is successful, tabbed text is returned. and we capture the fourth
# element as the requested mapped ID.
unlist(strsplit(content(response), "\\s+"))

# If the query can't be fulfilled because of a problem with the server, a
# WebPage is returned. But the server status is also returned and we can check
# the status code. I have lately gotten a few  "503" status codes: Server Not
# Available...

if (response$status_code == 200) { # 200: oK
    uniProtID <- unlist(strsplit(content(response), "\\s+"))[4]
    if (is.na(uniProtID)) {
        warning(paste("UniProt ID mapping service returned NA.",
                      "Check your RefSeqID."))
    }
} else {
    uniProtID <- NA
    warning(paste("No uniProt ID mapping available:",
                  "server returned status",
                  response$status_code))
}

# Let's see what we got...
uniProtID  # This should be "P39678" (or NA if the query failed)

# This is how "screenscraping" data from public Websites works in principle.


# ==== NCBI EUtils

# Next, we'll retrieve data from the various NCBI databases.

# It is unreasonably difficult to screenscrape the NCBI site directly, since the
# actual page contents are dynamically loaded via AJAX and not present in the
# page source that your browser makes available. This may be intentional, or
# just overengineering. NCBI offers a subset of their data via the EUtils API
# and that works well enough, but some of the data that is available on the Web
# browser pages is not served to a program.

# NCBI's eutils API returns data in XML format. I was able to simplify many,
# many lines of very techical legacy code from this material by using a new,
# convenient interface to EUtils that has recently been published on CRAN.
# Thanks to the author, Gerhard Schöfl, for the reutils package:

if (!require(reutils, quietly=TRUE)) {
    install.packages("reutils")
    library(reutils)
}

# Have a look at the following URL in your browser to see what e EUtils response
# looks like in principle:

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=NP_010227

# Fetch the UID associated with a RefSeqID
( UID <- esearch("NP_010227", "protein") )

# post the UID to the entrez History server and obtain a key
( eHist <- epost(UID, "protein") )

# Fetch this as an R data frame (mode "parsed"):
x <- content(esummary( eHist, version = "1.0"), "parsed")
str(x)
x$TaxId

( p <- efetch(eHist, retmode = "text", rettype = "fasta") )
# retrieve the content from the efetch object
content(p)
dbSanitizeSequence(content(p))

# Get the scientifc name for the tax-ID
( t <- efetch(x$TaxId, "taxonomy") )

# This is an XML document - we can access the node we want like so:
t$xmlValue("/TaxaSet/Taxon/ScientificName")

# So far so good. Since we need to do this a lot, we need to roll all of this
# into a function. I have added the function to the .utilities.R script. Study:
dbFetchProteinData

# Try:
x <- dbFetchProteinData("NP_010227")    # MBP1_SACCE
x <- dbFetchProteinData("XP_016612325") # MBP1_SPIPU
x <- dbFetchProteinData("NP_593032")    # MBP1_SCHPO

# Now, to put this to work we'll import all sequences that the PSI BLAST search
# has retrieved previously for our APSES domain protein seearches into our
# database. I have provided the file PSI-BLAST_hits.txt that contains the
# RefSeqIDs:
PSIBLASThits <- read.csv("PSI-BLAST_hits.txt",
                         header = FALSE,
                         stringsAsFactors = FALSE)
colnames(PSIBLASThits) <- "ID"
PSIBLASThits

# Now we can iterate over this list, fetch all the data, and add it to our
# datamodel. Let's be cautious though and not add this directly to our database,
# but to intermediate tables that we can sanity check before we merge them with
# myDB. Consider the code carefully, line by line.
#
# Note: the first row of the protein table gets an ID based on the IDs
# present in myDB, the subsequent ones get IDs based on the current contents of
# the intermediate table we are building.
#
# Note: the loop uses three control structures which you may not have
# encountered yet: conditional expressions using any(), and all(), and the
# next() statement which skips to the next iteration of a loop.

tmpDB <- dbInit()

for (refSeqID in PSIBLASThits$ID) {

    if (any(myDB$protein$RefSeqID == refSeqID)) { # we already got
                                                  # this RefSeqID
        next()
    }
    # get the data for one RefSeqID
    protData <- dbFetchProteinData(refSeqID)
    if (protData$version != myDB$version) {
        stop("Panic: datamodel versions do not match.")
    }
    # Create an ID for the protein table
    if (nrow(tmpDB$protein) == 0) {  # Create first ID
        protData$protein$ID <- dbCreateID(myDB$protein$ID)
    } else { # Create subsequent IDs
        protData$protein$ID <- dbCreateID(tmpDB$protein$ID)
    }
    # bind the new protein data to tmpDB
    tmpDB$protein <- rbind(tmpDB$protein, protData$protein)
    # bind the new taxonomy data to tmpDB
    existingTaxIDs <- unique(c(myDB$taxonomy$ID, tmpDB$taxonomy$ID))
    if (all(protData$taxonomy$ID != existingTaxIDs)) { # we don't have this
                                                       # taxonomy ID yet
        tmpDB$taxonomy <- rbind(tmpDB$taxonomy, protData$taxonomy)
    }
}

# Examine the new protein and taxonomy tables in tmpDB carefully:
# - are all the data elements present?
# - are the IDs in the protein table consecutive with the IDs in myDB?
# - are the IDs unique?

# If there is a problem, it needs to be fixed before you continue. If all is
# good, we can merge the two new tables into myDB.

myDB$protein <- rbind(myDB$protein, tmpDB$protein)
myDB$taxonomy <- rbind(myDB$taxonomy, tmpDB$taxonomy)

# If you have look at the taxonomy table carefully, you will have noticed that
# some organisms appear twice: once at the species level and once at the strain
# level. That is because I have originally annotated the reference proteins at
# the species level, the database records retrieved via RefSeqIDs are however
# annotated with the actual strain, from which the protein was sequenced. That's
# oK - as long as we keep it in mind when we select by species. In the "real"
# world, I would now go back and reannotate the reference database.


# Once you are sure the new data is correct, save the database into a new
# version.
save(myDB, file = "myDB.05.RData")



# ==============================================================================
#      PART NINE: COMPUTING OVER MSAs
# ==============================================================================

# ==  INFORMATION
# Shannon information is calculated as the difference between
# expected and observed entropy, where entropy is the negative sum over
# probabilities, times the log of those probabilities:

# In this section we compute Shannon information scores for aligned positions of
# the APSES domains we have collected with our PSI-BLAST search, and plot the
# values.

# === Step one: APSES domain MSA

# We start by collecting the APSES domain sequences from a pairwise alignment of
# the MBP1_SACCE domain with all others in the database. Note that this is NOT
# the best way to do such a domain alignment - using a HMM of annotated domains
# via SMART, Pfam, or ProSite, for example would be the Right Way - but it will
# work for our case (using only code that we already have), since the pairwise
# similarities are quite high. In a "final" state of our database, we would get
# the domains from their annotations, but all the new domains have not been
# annotated (yet).

# The reference domain for the alignment ...
aaRef_APSES <- AAString(dbGetFeatureSequence(myDB,
                                             "MBP1_SACCE",
                                             "APSES fold"))

APSESdomains <- character(nrow(myDB$protein))

for (i in 1:length(myDB$protein$ID)) {
    ali <-  pairwiseAlignment(
        aaRef_APSES,
        AAString(myDB$protein$sequence[i]),
        type = "overlap",
        substitutionMatrix = "BLOSUM62",
        gapOpening = 10,
        gapExtension = 0.5)
    # get sequence, but drop hyphens for indels
    APSESdomains[i] <- dbSanitizeSequence(as.character(ali@subject))
    names(APSESdomains)[i] <- myDB$protein$name[i]
}

# Now calculate the MSA
apsesSet <- AAStringSet(APSESdomains)
apsesSet@ranges@NAMES <- names(APSESdomains)
apsesMSA <- msaClustalOmega(apsesSet, order = "aligned")

# confirm
writeSeqSet(apsesMSA, format = "ali")

# === Step two: accessing the SeqSet column by column:

tmp <- as.character(apsesMSA)
apsesMat <- matrix(character(nchar(tmp[1]) * length(tmp)), nrow = length(tmp))
for (i in 1:length(tmp)) {
    apsesMat[i, ] <- unlist(strsplit(tmp[i], ""))
}
row.names(apsesMat) <- names(tmp)
colnames(apsesMat) <- 1:ncol(apsesMat)
# apsesMat contains the aligned characters, one character per element and we can
# easily extract them by row or column.
head(apsesMat)

apsesMat[ , 67]
table(apsesMat[ , 67])


# === Step three: define a function to calculate entropy

# Study the function code, then execute it to define it.
Htable <- function(tab) {
    # calculate shannon entropy for a frequency table tab.
    # Value: entropy H in log_2 units
    #
    # Note: we are ignoring gap characters in tab, since they do not
    # correspond to observations.
    #
    # Note: we are not correcting for small sample sizes
    # here. Thus if there are a large number of gaps in
    # the alignment, this will look like small entropy
    # since only a few amino acids are present. In the
    # extreme case: if a position is only present in
    # one sequence, that one amino acid will be treated
    # as 100% conserved - zero entropy. Sampling error
    # corrections are discussed eg. in Schneider et al.
    # (1986) JMB 188:414

    tab <- tab[names(tab) != "-"]
    tab <- tab/sum(tab)
    tab <- tab * (log(tab)/log(2))
    return(-sum(tab)) # return Shannon entropy
}

# Try:
Htable(table(apsesMat[ , 67]))

# _information_ is defined as the difference in entropies between the expected
# character frequencies and the the observed character frequencies (from an
# aligned column). Expected frequencies can be taken as 1/20 (all amino acids
# equally likely) ...

fAAuniform <- rep(1/20, 20)
names(fAAuniform) <- unlist(strsplit("ARNDCQEGHILKMFPSTWYV", ""))

( Huniform <- Htable(fAAuniform) ) # 4.3219 - this is the maximum
                                   # entropy for 20 characters


# Another useful set of expected frequencies can be taken from the UniProt
# database frequency statistics. cf. http://www.ebi.ac.uk/uniprot/TrEMBLstats -
# Section 5

fAAuniprot <- c(
    "A" =  8.92,
    "R" =  5.66,
    "N" =  3.95,
    "D" =  5.42,
    "C" =  1.26,
    "Q" =  3.84,
    "E" =  6.15,
    "G" =  7.16,
    "H" =  2.22,
    "I" =  5.70,
    "L" =  9.83,
    "K" =  5.06,
    "M" =  2.40,
    "F" =  3.93,
    "P" =  4.87,
    "S" =  6.81,
    "T" =  5.59,
    "W" =  1.30,
    "Y" =  2.95,
    "V" =  6.81
)

( Huniprot <- Htable(fAAuniprot) ) # 4.163

# What is the information in a column of alignment?
Huniform - Htable(table(apsesMat[ , 67])) # in bits


# With his, we can calculate information for each position of our multiple
# sequence alignment. Note that this code does not implement any sampling bias
# correction, so positions with a large numbers of gaps will receive
# artificially high scores due to small sample bias. In the extreme case, a
# single observed character is deemed to be 100% conserved, has an entropy of 0,
# and thus an information of 4.322:
Huniform - Htable(table("A"))

# ... and this will happen whenever we have mostly gaps:
Huniform - Htable(table(c("A", "-", "-", "-", "-", "-", "-")))

# result vector of Information values
vI <- numeric(ncol(apsesMat))

for (i in 1:ncol(apsesMat)) {
    vI[i] = Huniform - Htable(table(apsesMat[ , i]))
}

plot(vI)   # Let's improve this view a bit by adding color values ...
# evaluate
quantile(vI)
hist(vI)

# you can see that we have quite a large number of columns with the same,
# high value ... what are these?

which(vI > 4)
apsesMat[ , which(vI > 4.2)]

# And what is in the columns with low values?
apsesMat[ , which(vI < 1.1)]

# plot with colored bars
IP <- (vI-min(vI))/(max(vI) - min(vI) + 0.0001) # normalize range
nCol <- 15
IP <- floor(IP * nCol) + 1
spect <- colorRampPalette(c("#DD0033", "#00BB66", "#3300DD"), bias=0.6)(nCol)
Icol <- vector()
for (i in 1:length(vI)) {
    Icol[i] <- spect[ IP[i] ]
}

plot(1, 1,
     xlim=c(0, ncol(apsesMat)), ylim=c(-0.5, 5),
     type="n", bty="n",
     xlab="position in alignment", ylab="Information (bits)")

# draw values as colored rectangles: height is information and color is coded to
# information in the vector Icol
for (i in 1:ncol(apsesMat)) {
    rect(i, 0, i+1, vI[i], border=NA, col=Icol[i])
}

# As you can see, some of the columns reach very high values, but they are not
# contiguous in sequence. Most of these "conserved" columns are gaps in the
# alignment, but some correspond to actually very highly conserved rsidues.
# Overall, such plots of values along sequence coordinates are not all that
# useful. It becomes _a lot_ (!) more interesting, when we map these values to
# 3D structure. This will be our task in the next unit.

# Before you close this session, save apsesMat and vI for the Structure.R unit:
save(apsesMat, vI, file = "infoValues.RData")



# This concludes the R code for the Sequence unit.
#  - myDB should be easily recreatable from version 0.5, which includes
#    all APSES domain proteins we found via PSI BLAST
#  - a number of tasks have been mentioned along the way, make sure
#    you have completed them.
#  - Work through the Exercises now.


# ==============================================================================
#      PART TEN: EXERCISES
# ==============================================================================
# Please work through the following exercise questions.
# I may ask you to hand in your solutions for credit at
# some point in the course.

# SEQUENCE-1:
# Write a regular expression that you can use with grep() to determine if a
# string contains both letters and digits. Apply this to find "names" in the
# protein table that do _not_ have letters and digits.

# SEQUENCE-2:
# aaindex[[459]] lists the amino acid frequencies for cytoplasmic proteins.
# aaindex[[460]] lists the frequencies for nuclear proteins. Transcription
# factors are nuclear proteins. (a) evaluate the difference between these two
# tables (compute log ratios): which amino acids are realtively more/less
# abundant in the nucleus? Does this make sense according to your understanding
# of biology? (b) redo the enrichment plot with nuclear background frequencies.
# What does the plot now mean? Does this make sense?

# SEQUENCE-3:
# Produce a dotplot of the ankyrin domain region of MBP1_SACCE against itself.
# Interpret.

# SEQUENCE-4:
# Rewrite the code of plotProtein() so that the N-termini of annotated APSES-,
# Ankyrin-, and Swi6 domains are aligned. Hint: Develop this first "manually",
# with a copy of MBP1_SCHPO (that seems to be the shortest protein). From all
# squences to be plotted, first find the maximum coordinate for the respective
# domain that needs to be accommodated, then move the relevant annotation start-
# and endpoints. Then overwrite the "inserted" sequence with a white box as the
# last "annotation" to be plotted.

# SEQUENCE-5:
# Calculate an MSA for MBP1_... proteins' Swi6 domains. Write the sequence set
# to an MFA file and compare the MAFFT and TCoffee alignments at the EBI.

# SEQUENCE-6:
# Calculate an MSA for MBP1_... proteins' Ankyrin domains.
# - Start the sequence range to be aligned from the smallest Start coordinate
#   for the feature in a sequence. Hint: use the min() function.
# - End the sequence range at the largest End coordinate. Hint: use the max()
#   function. Write the sequence set to an MFA file and compare the MAFFT and
#   TCoffee alignments at the EBI.

# SEQUENCE-7:
# I am not happy with our code to infer protein "names" from FASTA headers in
# the $Title attribute of the NCBI entrez response. Some proteins are named
# "protein", "factor", or "partial". How would you improve this code? Can the
# regex you have devloped in the SEQUENCE-1 exercise help? Share your ideas
# on the mailing ist.

# SEQUENCE-8:
# Write an expression that subsets apsesMat to contain only those columns that
# are not "-" characters in a given row - e.g. the row that has rownames()
# "MBP1_SACCE"


# You are welcome to eMail if you need additional hints or help. I will respond
# if your mail conforms to the standards laid out in the two documents in the
# assigment Wiki Page footer. You may assume that I have refDB and myDB
# available, i.e. you do not need to include instructions on how to recreate
# a database in your MWE.



# [END]
