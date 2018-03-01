package alignmentprofile

import (
	a "github.com/hivdb/nucamino/types/amino"
)

type Gene string
type PositionalIndelScores map[int]([2]int)
type GenePositionalIndelScores map[Gene]PositionalIndelScores
type ReferenceSeqs map[Gene][]a.AminoAcid

// This stores the all the information needed to align a sequence
// against a particular reference: alignment parameters, positional
// in-del scores, and the reference sequences to be aligned against.
type AlignmentProfile struct {
	StopCodonPenalty         int
	GapOpeningPenalty        int
	GapExtensionPenalty      int
	IndelCodonOpeningBonus   int
	IndelCodonExtensionBonus int
	GeneIndelScores          GenePositionalIndelScores
	ReferenceSequences       ReferenceSeqs
}

func (profile *AlignmentProfile) SupportsPositionalIndelScores() bool {
	return profile.GeneIndelScores != nil
}

func (profile *AlignmentProfile) SupportsPositionalIndelScoresFor(gene Gene) bool {
	if !profile.SupportsPositionalIndelScores() {
		return false
	}
	_, hasKey := profile.GeneIndelScores[gene]
	return hasKey
}

func (profile *AlignmentProfile) PositionalIndelScoresFor(g Gene) (PositionalIndelScores, bool) {
	if !profile.SupportsPositionalIndelScores() {
		return nil, false
	}
	scores, found := profile.GeneIndelScores[g]
	return scores, found
}

// An array of all the genes supported by this alignment profile.
func (profile *AlignmentProfile) Genes() []Gene {
	var genes = make([]Gene, 0, len(profile.ReferenceSequences))
	for gene, _ := range profile.ReferenceSequences {
		genes = append(genes, gene)
	}
	return genes
}
