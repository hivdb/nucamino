// This package contains data structures that hold all the information
// required by the alignment algorithm, and procedures for serializing
// and deserializing that information.
package alignmentprofile

import (
	"fmt"
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

func (profile AlignmentProfile) asRaw() rawAlignmentProfile {
	var raw rawAlignmentProfile
	raw.StopCodonPenalty = profile.StopCodonPenalty
	raw.GapOpeningPenalty = profile.GapOpeningPenalty
	raw.GapExtensionPenalty = profile.GapExtensionPenalty
	raw.IndelCodonOpeningBonus = profile.IndelCodonOpeningBonus
	raw.IndelCodonExtensionBonus = profile.IndelCodonExtensionBonus

	raw.ReferenceSequences = make(map[string]string)
	for gene, aaSeq := range profile.ReferenceSequences {
		raw.ReferenceSequences[string(gene)] = a.WriteString(aaSeq)
	}

	if profile.GeneIndelScores != nil {
		raw.RawIndelScores = profile.rawIndelScores()
	}
	return raw
}

// Retrieve the positional indel scores for a Gene.
func (profile *AlignmentProfile) PositionalIndelScoresFor(g Gene) (PositionalIndelScores, bool) {
	scores, found := profile.GeneIndelScores[g]
	return scores, found
}


// Check that the profile isn't empty
func (profile AlignmentProfile) validate() error {
	if len(profile.ReferenceSequences) == 0 {
		return fmt.Errorf("Missing key: ReferenceSequence")
	}
	return nil
}
