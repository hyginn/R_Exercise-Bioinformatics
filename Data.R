# Data.R
#
# Purpose:  Introduction to Bioinformatics with R - Data
#
# Version: 1.0
#
# Date:    2016  12  30
# Author:  Boris Steipe (boris.steipe@utoronto.ca)
#
# V 1.0    First final version
#
# TODO:
#
#
# ==============================================================================

# ==============================================================================
#      PART ONE: SETUP
# ==============================================================================

# We will explore transcription factors in the newly sequenced, archaic fungus
# Spizellomyces punctatus, but we will start out from proteins that have been
# characterized in Saccharomyces cerevisiae, and nine other "reference" species
# that are well-distributed over the phylogenetic tree.

# Since we often need a shorthand notation of our species of interest to label
# sequences or other objects, I have included a function: bioCode() to extract
# letters from a binomial species name. Look at the function (typing a function
# name _without_ parentheses prints the function code):

biCode

# Then use it to display the "biCode" of "Spizellomyces punctatus":
biCode("Spizellomyces punctatus")

# I will usually use these bicode shorthand labels to refer to species.



# ==============================================================================
#      PART TWO: A PROTEIN DATABASE IMPLEMENTED AS A LIST OF DATAFRAMES
# ==============================================================================

# Every entity in our datmodel can be implemented by a data frame. I will refer
# to these data frames as "tables" as they are generally called this in
# relational databases. Thus "table" is not a special data structure of R, but a
# name we give to a data frame that we use to hold entity data in a database.
#
# To keep things organized, we create a list, that contains the entity tables.
# Here is a definition for a table for protein data as described in the data
# model schema:
#
db <- list()      # we'll call our database "db" and start with an empty list

# We create a data frame that implements the "protein" entity, with a single
# row of data, and add it to the list. Note that we create the table directly in
# the list, simply by referencing its (as yet non-existent) name after the "$"
# operator.
db$protein <- data.frame(ID = "my_pro_1",
                         name = "Mbp1",
                         RefSeqID = "NP_010227",
                         UniProtID = "P39678",
                         taxonomy.ID = as.integer(4932),
                         sequence = "...",              # just a placeholder
                         stringsAsFactors = FALSE)

objectInfo(db)

# Next, we create the taxonomy table and add that to the list. Note that as
# before, we simply assign a value (a data frame with two columns, called "ID",
# and "species") to a new object (...$taxonomy) of our list (db).

db$taxonomy <- data.frame(ID = 4932L,
                          species = "Saccharomyces cerevisiae",
                          stringsAsFactors = FALSE)

objectInfo(db)

# We can address and select these items with the "$" or "[[]]" operators:
# ... a data frame:
db$taxonomy
db[["taxonomy"]]

# a column from that data frame:
db$taxonomy$species
db[["taxonomy"]]$species
db$taxonomy[["species"]]


# Let's add a second protein: protein data ...
myRow <- data.frame(ID = "my_pro_2",
                    name = "Res2",
                    RefSeqID = "NP_593032",
                    UniProtID = "P41412",
                    taxonomy.ID = as.integer(4896),
                    sequence = "...",              # again, just a placeholder
                    stringsAsFactors = FALSE)

# ... in the protein data frame, wich we we rbind() to the existing table in "db":
( db$protein <- rbind(db$protein, myRow) )


# Same thing for taxonomy:
myRow <- data.frame(ID = as.integer(4896),
                    species = "Schizosaccharomyces pombe",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

objectInfo(db)

# Here are some examples again of accessing the contents:
db$protein
db$protein$RefSeqID
db$protein[, "name"]
db$taxonomy
db$taxonomy$species
biCode(db$taxonomy$species)


# ==== Selecting by condition
#
# Let's look up information in one table, based on a condition in another table:
# what is the species name for the protein whose name is "Mbp1"?

# First, we need to retrieve the taxonomy.ID of the Mbp1 protein. This is the
# "foreign key" we need to retrieve the correct information in the taxonomy
# table. Doing this, we use the relationship between a key, which is attribute
# of one item, and the contents of another table. We find the key in a cell in
# the table: db$protein[<row>, <column>], where the element in <row> is that one
# which the value in the "name" column is Mbp1:

db$protein$name == "Mbp1" # TRUE FALSE

# The <column> is called "taxonomy.ID". Simply
# insert these two expressions in the square
# brackets.

db$protein[db$protein$name == "Mbp1", "taxonomy.ID"]

# Assign the taxonomy.ID value ...
( taxID <- db$protein[db$protein$name == "Mbp1", "taxonomy.ID"] )

# ... and fetch the species_name value from db$taxonomy

db$taxonomy[db$taxonomy$ID == taxID, "species"]


# ==== Preparing and entering sequence

# The Mbp1 sequence as copied from the NCBI Website
# (https://www.ncbi.nlm.nih.gov/protein/NP_010227) contains characters that are
# not part of the sequence: blanks, and sequence numbers. We need to pre-process
# such data. In general: all data that goes into a databse should run through
# checks that ensure consistency and integrity. In this case, we want to remove
# all characters that are not part of the sequence itself, including whitespace,
# and we want to change the one-letter code to uppercase. First, assign the data
# as copy/pasted from the NCBI webpage.

( mySeq <- "
1 msnqiysary sgvdvyefih stgsimkrkk ddwvnathil kaanfakakr trilekevlk
61 ethekvqggf gkyqgtwvpl niakqlaekf svydqlkplf dftqtdgsas pppapkhhha
121 skvdrkkair sastsaimet krnnkkaeen qfqsskilgn ptaaprkrgr pvgstrgsrr
181 klgvnlqrsq sdmgfprpai pnssisttql psirstmgpq sptlgileee rhdsrqqqpq
241 qnnsaqfkei dledglssdv epsqqlqqvf nqntgfvpqq qssliqtqqt esmatsvsss
301 pslptspgdf adsnpfeerf pgggtspiis miprypvtsr pqtsdindkv nkylsklvdy
361 fisnemksnk slpqvllhpp phsapyidap idpelhtafh wacsmgnlpi aealyeagts
421 irstnsqgqt plmrsslfhn sytrrtfpri fqllhetvfd idsqsqtvih hivkrksttp
481 savyyldvvl skikdfspqy rielllntqd kngdtalhia skngdvvffn tlvkmgaltt
541 isnkegltan eimnqqyeqm miqngtnqhv nssntdlnih vntnnietkn dvnsmvimsp
601 vspsdyityp sqiatnisrn ipnvvnsmkq masiyndlhe qhdneikslq ktlksisktk
661 iqvslktlev lkesskdeng eaqtnddfei lsrlqeqntk klrkrliryk rlikqkleyr
721 qtvllnklie detqattnnt vekdnntler lelaqeltml qlqrknklss lvkkfednak
781 ihkyrriire gtemnieevd ssldvilqtl iannnknkga eqiitisnan sha
//
" ) # "\n" means: line-break

# To remove all characters that are not
# letters we use gsub(). gsub() "globally substitutes" a string for any stringThe "substitution" we specify is: choose every character NOT in the
# range between a-z or between A-Z and replace it with "" - i.e. the empty
# string. The syntax of "[^a-zA-Z]" is the standard "regular expression syntax:
# the square brackets in _this_ context hold values for a "character class" (any
# character to which the contents of the square brakets applies) and the caret -
# "^", inverts the definition, i.e. is a "not" operator.

( mySeq <- gsub("[^a-zA-Z]", "", mySeq) ) # replace all non-letters with ""

# Now to change the sequence to uppercase. R has the functions toupper()
# and tolower().

toupper(mySeq)

# I have combined these two steps into a function in the .utilities.R script.
# View:

dbSanitizeSequence

# (Remember: executing a function name without parentheses displays the
#  function code.)

# Try the function:

dbSanitizeSequence("f123  a$^&&*)m. i	l马马虎虎 é yßv+w")
dbSanitizeSequence(mySeq)


# === Updating the database ====================================================

# Modifying a field

# Here is code that modifies the sequence field in the protein table of the
# database:
#

str(db$protein) # before

( select <- db$protein$name == "Mbp1" )
( value <- dbSanitizeSequence(mySeq) )
db$protein$sequence[select] <- value

str(db$protein) # after

# Analyze the expression. Note how we specify an element of a vector (column) in
# a data frame in a list using a logical expression. And note how we assign the
# output (return value) of a function. As far as database coding goes this is
# pretty minimal - there is no error checking done at all. In particular: can we
# really guarantee that the name "Mbp1" is unique in the protein table? No! We
# never required it to be unique. This is a check we need to perform so
# frequently that we will encapsulate it in a function. Again, the function
# implements both sanity- and value checks.

dbConfirmUnique

# try this:
dbConfirmUnique(c("TRUE", "FALSE")) # Nope" "TRUE" is not TRUE
dbConfirmUnique(c(TRUE, FALSE))     # Yes
dbConfirmUnique(c(TRUE, TRUE))      # No: more than one
dbConfirmUnique(c(FALSE, FALSE))    # Also not good: less than one

# ==== Continuing to update the sequence field.

# Here is the update to the sequence field of Res2 but using our
# confirmUnique() function

mySeq <- "
1 maprssavhv avysgvevye cfikgvsvmr rrrdswlnat qilkvadfdk pqrtrvlerq
61 vqigahekvq ggygkyqgtw vpfqrgvdla tkykvdgims pilsldideg kaiapkkkqt
121 kqkkpsvrgr rgrkpsslss stlhsvnekq pnssisptie ssmnkvnlpg aeeqvsatpl
181 paspnallsp ndntikpvee lgmleapldk yeeslldffl hpeegripsf lyspppdfqv
241 nsvidddght slhwacsmgh iemiklllra nadigvcnrl sqtplmrsvi ftnnydcqtf
301 gqvlellqst iyavdtngqs ifhhivqsts tpskvaaaky yldcilekli siqpfenvvr
361 lvnlqdsngd tslliaarng amdcvnslls ynanpsipnr qrrtaseyll eadkkphsll
421 qsnsnashsa fsfsgispai ispscsshaf vkaipsissk fsqlaeeyes qlrekeedli
481 ranrlkqdtl neisrtyqel tflqknnpty sqsmenlire aqetyqqlsk rlliwlearq
541 ifdlerslkp htslsisfps dflkkedgls lnndfkkpac nnvtnsdeye qlinkltslq
601 asrkkdtlyi rklyeelgid dtvnsyrrli amscginped lsleildave ealtrek
"

( select <- dbConfirmUnique(db$protein$name == "Res2") )
( value <- dbSanitizeSequence(mySeq) )
db$protein$sequence[select] <- value

objectInfo(db$protein)

# === Adding a record

# Adding a record is easy in principle, simply defining the values we get from
# the NCBI or EBI databases ... except for the ID field. That is a field we need
# to define internally, and we absolutely need to make sure this is done right.
# Getting the primary- or foreign keys wrong or inconsistent is a recipe for
# disaster in a database. Once again, we'll write a small function: to define an
# ID for a new entry to a database table. The structure of an ID in our database
# will be a little more complex than just making a unique integer. Database
# purists may say that IDs should not contain _any_ semantics, because we should
# never have to interpret them in a way other than merely being a reference to
# an entry in a database table. But that only holds true if all the work of
# selecting and cross-referencing and sanity checking is done by the database
# engine, e.g. in a mySQL relational database. In our case, we are basically
# doing things by hand to learn how such a database works in principle. Adding a
# bit of sementics to our IDs is thus the Right Way. There are two elements that
# we will add to IDs:

#    IDs will have a "namespace" so we can distinguish who created an entry: if
#    this was me, we will use "ref" (for reference). If it was you, the
#    namespace prefix will be "my".

#    IDs will also contain a "tableCode", that specifies in which table they
#    were used. If we need to troubleshoot an error, it may be helpful to
#    understand which table an ID came from.

dbCreateID

# Let's analyze the first part of the code: creating a tableCode. missing() is a
# function we often use inside functions. It checks whether a parameter value
# has been defined for an argument. In our case, the first if() { ... } block
# checks for tableCode. If this is not missing(), then we continue right away:
# the user has explicitly requested a specific code. But if a value for
# tableCode is missing(), we need to create a value. As explained in the code,
# we have three options: either we get the code from the first element, if there
# is no first element yet we try to get the code from an attribute of the
# column, and if there is none defined, we assign a default value of "xxx".

# The second part of the function creates the actual ID. The
# first expression requires a bit of care to analyze. At its centre is the
# expression:
#    grep(sprintf("^%s_", ns), x)

# This follows the pattern: grep(<pattern>, <object>) - grep returns the index
# of the elements that contain a particular pattern. This pattern is constructed
# from a sprintf() expression which uses "ns": for example if the string is
# "ref", the pattern turns out as "^ref_". In this way, we can search for
# strings "dynamically - i.e. we don't need to know at the time we write the
# code what pattern we might be looking for. But what does "^ref_" mean? The
# "ref_" part is simply the literal string "ref_". But the caret in _this_
# context of a regular expression does not mean "NOT", but it means:
# "at-the-beginning-of-the-string". You need to read regular expressions
# carefully - the meaning of characters inside and ourside of square brackets
# may not be the same. For example "$" is the literal "$" character inside
# brackets, but outside of brackets it is the opposite to "^" and means:
# "at-the-end-of-the-string". Consider:

test <- c("ref_pref_1", "my_pref_1", "my_pref_2", "ref_pref_2")
ns <- "ref"

grep("ref", test)
sprintf("%s_", ns)
grep(sprintf("%s_", ns), test)
sprintf("^%s_", ns)
grep(sprintf("^%s_", ns), test)
test[grep(sprintf("^%s_", ns), test)]
test[-(grep(sprintf("^%s_", ns), test))]

# Ok. But what's the point? This bit of code separates out the target namespace
# from all IDs that the column might already contain. Through this we maintain
# the namespaces in the column completely independent of each other when we
# create new IDs. The variable "iIDs" (index of IDs) contains only indices for
# the IDs we actually need to consider when incrementing the <integer> part of
# the ID by 1 in the last expression. We separate out the numeric part from all
# these keys by removing the prefix, then look for the largest number among
# these with max(), and use an integer that is one larger than that to create a
# guaranteed unique key.

# After studying the function code, predict the results:
dbCreateID(character())
dbCreateID(character(), tableCode = "tbl")
dbCreateID(character(), ns = "tmp")
dbCreateID(sprintf("my_tbl_%d", 1:4))
dbCreateID(sprintf("ref_tbl_%d", 1:5))
dbCreateID(sprintf("ref_tbl_%d", 1:8), ns = "new")


# Using this function, here is the complete code to add one new item into the
# protein table:

mySeq <- "
1 msgdktifka tysgvpvyec iinnvavmrr rsddwlnatq ilkvvgldkp qrtrvlerei
61 qkgihekvqg gygkyqgtwi pldvaielae ryniqgllqp itsyvpsaad spppapkhti
121 stsnrskkii padpgalgrs rratsietes evigaapnnv segsmspsps dissssrtps
181 plpadrahpl hanhalagyn grdannhary adiildyfvt enttvpslli npppdfnpdm
241 sidddehtal hwacamgrir vvklllsaga difrvnsnqq talmratmfs nnydlrkfpe
301 lfellhrsil nidrndrtvf hhvvdlalsr gkphaaryym etminrlady gdqladilnf
361 qddegetplt maararskrl vrlllehgad pkirnkegkn aedyiieder frsspsrtgp
421 agielgadgl pvlptsslht seagqrtagr avtlmsnllh sladsydsei ntaekkltqa
481 hgllkqiqte iedsakvaea lhheaqgvde erkrvdslql alkhainkra rddlerrwse
541 gkqaikrarl qaglepgals tsnatnapat gdqkskddak sliealpagt nvktaiaelr
601 kqlsqvqank telvdkfvar areqgtgrtm aayrrliaag cggiapdevd avvgvlcell
661 qeshtgarag aggerddrar dvammlkgag aaalaanaga p
"

myRow <- data.frame(ID = dbCreateID(db$protein$ID),
                    name = "UMAG_1122",
                    RefSeqID = "XP_011392621",
                    UniProtID = "A0A0D1DP35",
                    taxonomy.ID = as.integer(5270),
                    sequence = dbSanitizeSequence(mySeq),
                    stringsAsFactors = FALSE)
db$protein <- rbind(db$protein, myRow)

myRow <- data.frame(ID = as.integer(5270),
                    species = "Ustilago maydis",
                    stringsAsFactors = FALSE)
db$taxonomy <- rbind(db$taxonomy, myRow)

objectInfo(db)

# Should we also put this all together into a function? I have not done this
# because copying this actual code and changing only the parameter values is not
# more involved than putting the values into the argument list of a function.
# For a "real" database application, we would ceratinly have a function to add
# values - but it would do much more than simply adding values into tables: it
# would implement sanity and consistency checks for all items. This is beyond
# the scope of this simple example, but of course you could try your hand at
# writing such a dbAddProtein() function.

# Adding items with such raw code templates without error checking can certainly
# lead to mistakes. If we notice an error, often we can just overwrite the
# offending field with correct data. But sometimes it will be easier (and more
# robust) to delete the erroneous entry and add the correct one. For example, if
# your code is in a script, and you realize the entry had an error, I would not
# "patch" the error in the script by changing a value in the database, but I
# would delete the entry by enetering the appropriate command into the Console,
# correct the script, and execute the corrected block. That way the mistake is
# fixed at its source and the script become the authoritative source of correct
# commands.

# Removing a row from a datframe is trivial: just overwrite the dataframe with
# a selection statement in which the unique selection of the offending row is
# inverted:

# create an erroneous entry
myRow <- data.frame(ID = dbCreateID(db$protein$ID),
                    name = "nonesuch",
                    RefSeqID = "NP_000000",
                    UniProtID = "A123456",
                    taxonomy.ID = as.integer(999999),
                    sequence = dbSanitizeSequence("very sequense much wow",
                                                  strictAA = FALSE),
                    stringsAsFactors = FALSE)
( db$protein <- rbind(db$protein, myRow) )

# Make a logical vector that identifies it
( select <- dbConfirmUnique(db$protein$name == "nonesuch") )
!select    # its logical inversion

str(db)    # before
db$protein <- db$protein[ ! select, ]  # overwrite the table with a copy
                                       # without the selected record
str(db)    # after

# Note: if you delete records "by hand" you need to be careful that you do
# not remove keys that are used as foreign keys in another table - if there
# are such dependencies, you need to update the other table(s) too. "Real"
# database systems include such dependencies in the creation instructions
# of the table schema: "on delete cascade ..."


# ==============================================================================
#      PART THREE: THE "REFERENCE DATABASE"
# ==============================================================================

# I have implemented the protein and taxonomy tables in a "reference database",
# added a few other tables that we will encounter later, and populated the
# database with MBP1 orthologues from ten fungal species. The database was
# created through the createRefDB.R script, which is included in the project
# files.

objectInfo(refDB)



# ==============================================================================
#      PART FOUR: EXERCISES
# ==============================================================================
# Please work through the following exercise questions.
# I may ask you to hand in your solutions for credit at
# some point in the course.



# DATA-1:
# Write an expression that returns the names of the objects in refDB



# DATA-2:
# Write a for-loop that lists the names of all objects that each of the objects
# in refDB contains. (e.g. the table "taxonomy" in refDB contains the objects
# "ID" and "species")



# DATA-3:
# Now enhance that loop to print the name of the object, the number of objects
# it contains, and their names. Here is an example for one line of the output:
#             taxonomy	(2)	ID	species



# DATA-4:
# Write a function that returns the name and species for all proteins in a
# database that is structured like refDB, whose name matches a regular
# expression pattern that you can pass as a parameter, for  each species whose
# name contains a pattern that you also pass as parameter. Make sure that no
# error is thrown if either parameter does not result in a match. Return the
# empty string "" if no match is found.

# Hint: develop this step by step.
# Hint: use grep() or grepl() as appropriate
# Hint: use the operator %in% to retrieve matches against more than one value
#       example:
#          a <- c(3, 7, 1)
#          b <- 1:10
#          b %in% a     # return TRUE for any b that is also in a
#          b[b %in% a]

#  Here is an example function call and output.
#
# getProteinByRegex(DB = refDB, name = "[YU]", species = "cc")
#
#         name                 species
# 7 MBP1_CRYNE Cryptococcus neoformans
# 8 MBP1_PUCGR       Puccinia Graminis
#



# DATA-5:
# Change your function so that if name or species are not specified, all
# names or species, respectively, will be matched by default (in case the
# function does not behave in this way already).



# You are welcome to eMail if you need additional hints or help. I will respond
# if your mail conforms to the standards laid out in the two documents in the
# assigment Wiki Page footer. You may assume that I have refDB available, i.e.
# you do not need to include instructions on how to recreate refDB in your MWE.


# [END]
