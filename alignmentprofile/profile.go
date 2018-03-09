package alignmentprofile

import (
	a "github.com/hivdb/nucamino/types/amino"
	"sort"
)

type Gene string
type PositionalIndelScores map[int]([2]int)
type GenePositionalIndelScores map[Gene]PositionalIndelScores
type ReferenceSeqs map[Gene][]a.AminoAcid

// This stores the all the information needed to align a sequence to a
// reference: reference sequences, alignment parameters, and
// positional indel scores.
type AlignmentProfile struct {
	StopCodonPenalty         int
	GapOpeningPenalty        int
	GapExtensionPenalty      int
	IndelCodonOpeningBonus   int
	IndelCodonExtensionBonus int
	GeneIndelScores          GenePositionalIndelScores
	ReferenceSequences       ReferenceSeqs
}

// An array of all the genes supported by this alignment profile.
func (profile AlignmentProfile) Genes() []Gene {
	var genes = make([]Gene, 0, len(profile.ReferenceSequences))
	for gene, _ := range profile.ReferenceSequences {
		genes = append(genes, gene)
	}
	return genes
}

func (profile AlignmentProfile) rawIndelScores() map[string][]rawIndelScore {
	result := make(map[string][]rawIndelScore)
	for gene, positionalIndelScores := range profile.GeneIndelScores {
		rawIndelScores := make([]rawIndelScore, 0, len(positionalIndelScores))
		for posKey, scores := range positionalIndelScores {
			var pos int
			var kind string
			// Interpret the position key in the positional indel
			// scores map: the magnitude is the position, negative
			// indicates deletion.
			if posKey >= 0 {
				pos = posKey
				kind = "ins"
			} else {
				pos = -posKey
				kind = "del"
			}
			rawScore := rawIndelScore{
				Kind:     kind,
				Position: pos,
				Open:     scores[0],
				Extend:   scores[1],
			}
			rawIndelScores = append(rawIndelScores, rawScore)
		}
		sort.Sort(byPositionAndKind(rawIndelScores))
		result[string(gene)] = rawIndelScores
	}
	return result
}
