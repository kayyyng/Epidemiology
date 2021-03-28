setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
aaAnnotation = function(dnaSet, refseq, geneseq, pos, WTaa = NULL, searchLen = 20) {
  # Single amino acid variant annotation
  # Annotate single amino acid variants of a set of sequences aligned to refseq
  ## Input:
    ## dnaSet: the sequence set of interest aligned to a reference sequence (DNAStringSet)
    ## refseq: the reference sequence after multiple sequence alignment, may contain gaps (DNAStringSet)
    ## geneseq: the nucleotide sequence of the gene where the variation locates (DNAStringSet)
    ## pos: the amino acid residue position of the variant (integer)
    ## WTaa: (optional) the amino acid residue in the reference sequence, for assertion (character)
    ## searchLen: the number of base pair used for pairwise-alignment of geneseq to refseq (integer)
  ## Output: a dataframe with
    ## sequence: sequence ID
    ## (pos): Amino acid variant, "X" if ambiguous
  library(Biostrings)
  initPos = pos*3-2
  queryseq <- geneseq[[1]][initPos:(initPos+searchLen-1)]
  print(queryseq) # Observe the query sequence, adjust the length and position if gaps exist
  
  # Align query sequence to reference sequence
  x = pairwiseAlignment(queryseq, refseq, type = "global-local")
  aaInitPos = start(subject(x))
  
  # Assertion
  if (WTaa != "X") {
    codon = refseq[[1]][aaInitPos:(aaInitPos+2)]
    aa = translate(codon, no.init.codon = TRUE, )
    if (aa != AAString(WTaa)){
      print("Error: unexpected amino acid")
      print(paste("Expected ", WTaa, ", obtained ", aa, sep = ""))
      break
    }
  }
  
  # Translate the DNA set
  codonSet = subseq(dnaSet, start = aaInitPos, end = aaInitPos+2)
  aaSet = translate(codonSet, no.init.codon = TRUE, if.fuzzy.codon = "X")
  df = data.frame("sequence" = names(aaSet), aaSet)
  names(df)[2] = paste(WTaa, toString(pos), sep = "")
  write.csv(df, file = paste(toString(pos), "aaAnnotation.csv", sep=""), row.names = FALSE)
  return(df)
}

# Implement the function
dnaSet = readDNAStringSet("large_with_Australia_aligned_to_refseq.fasta")
refseq = readDNAStringSet("refseq_aligned.fasta")
geneseq = readDNAStringSet("S.fna")
df = aaAnnotation(dnaSet, refseq, geneseq, pos = 222, WTaa = "A")

