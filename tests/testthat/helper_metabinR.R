#Dummy fasta (3 sequences of 10bp length) temp file

fastaFile <- tempfile(fileext = ".fasta")
fastaStr <- paste0(
    ">S1\n",
    "ACGTACGTAA\n",
    ">S2\n",
    "ACGTACGTCC\n",
    ">S3\n",
    "AAAAAGGGGG\n"
)
write(fastaStr, file = fastaFile, sep = "")
