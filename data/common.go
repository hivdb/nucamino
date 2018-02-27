package data

import a "github.com/hivdb/nucamino/types/amino"

// A map of gene names to arrays of amino acids, for storing the
// reference sequences of genes of interest.
type SequenceMap map[string][]a.AminoAcid
