
    library("DNAcopy")
    library(devtools)
    library(CRISPRcleanR)



    # load library

    args <- commandArgs(TRUE)  # read custom library
    libraryfile <- args[1]
    rcfiles <- args[2]

    if (strcmp(libraryfile, "Yusa")){
        data(KY_Library_v1.0)
        librarydata = KY_Library_v1.0
    } else {
        librarydata = read.table(libraryfile)
    }


    # CRISPRcleanR function to apply to each file
    CRISPRcleanR <- function(x,ctrln) {
      fn = paste(x, sep="	")
      normANDfcs = ccr.NormfoldChanges(fn, min_reads = 30, libraryAnnotation = librarydata, ncontrols = ctrln, outdir = './') # libraryAnnotation = TKO
      gwSortedFCs = ccr.logFCs2chromPos(normANDfcs$logFCs, librarydata)
      correctedFCs = ccr.GWclean(gwSortedFCs, display = FALSE, label = x)
      x_correctedFC = correctedFCs$corrected_logFCs
      write.table(x_correctedFC, paste(sub(".readcount", "", x), "CCR" , sep="."), sep="	")
    }


    filelist <- read.table(rcfiles,stringsAsFactors = FALSE)

    # apply CRISPRcleanR function on each file in files
    for (row in 1:nrow(filelist)) {
        filename <- filelist[row, "FILE"]
        controln  <- filelist[row, "CTRLN"]
        CRISPRcleanR(filename,controln)
    }
    