package alignment

import (
	a "github.com/hivdb/nucamino/types/amino"
)

type Gene string

type PositionalIndelScores map[Gene]map[int]([2]int)

// This stores the all the information needed to align a sequence
// against a particular reference: alignment parameters, positional
// in-del scores, and the reference sequences to be aligned against.
type AlignmentProfile struct {
	StopCodonPenalty         int
	GapOpenPenalty           int
	IndelCodonOpeningBonus   int
	IndelCodonExtensionBonus int
	GeneIndelScores          PositionalIndelScores
	ReferenceSequences       map[Gene][]a.AminoAcid
}

func (profile *AlignmentProfile) SupportsPositionalIndelScores() bool {
	return profile.GeneIndelScores != nil
}

func (profile *AlignmentProfile) SupportsPositionalIndelScoresFor(gene Gene) bool {
	if ! profile.SupportsPositionalIndelScores() {
		return false
	}
	_, hasKey := profile.GeneIndelScores[gene]
	return hasKey
}

func (profile *AlignmentProfile) Genes() []Gene {
	var genes = make([]Gene, 0, len(profile.ReferenceSequences))
	for gene, _ := range profile.ReferenceSequences {
		genes = append(genes, gene)
	}
	return genes
}
