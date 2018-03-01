package alignment

import (
	a "github.com/hivdb/nucamino/types/amino"
	"testing"
)

var exampleProfile = AlignmentProfile{
	ReferenceSequences: map[Gene][]a.AminoAcid{
		"A": a.ReadString("T"),
		"B": a.ReadString("V"),
	},
}

var exampleProfileWithPositionalIndelScores = AlignmentProfile{
	ReferenceSequences: map[Gene][]a.AminoAcid{
		"A": a.ReadString("T"),
		"B": a.ReadString("V"),
	},
	GeneIndelScores: PositionalIndelScores{
		"A": map[int][2]int{
			0: [2]int{1, 2},
			1: [2]int{3, 4},
		},
	},
}

func TestGenes(t *testing.T) {
	genes := exampleProfile.Genes()
	expected := [2]Gene{"A", "B"}

	if len(genes) != len(expected) {
		t.Errorf("Expected %v to be %v", genes, expected)
	}

	for _, exp := range expected {
		found := false
		for _, g := range genes {
			if g == exp {
				found = true
			}
		}
		if !found {
			t.Errorf("Missing key: %v", exp)
		}
	}

}

func TestPositionalIndelScoreSupportChecking(t *testing.T) {
	if exampleProfile.SupportsPositionalIndelScores() {
		t.Errorf("This example profile doens't have positional indel scores")
	}
	if !exampleProfileWithPositionalIndelScores.SupportsPositionalIndelScores() {
		t.Errorf("This example profile does have positional indel scores")
	}
}

func TestPositionalIndelScoresByGene(t *testing.T) {
	p := exampleProfileWithPositionalIndelScores
	if !p.SupportsPositionalIndelScoresFor("A") {
		t.Errorf("Expected this profile to have scores for 'A'")
	}
	if p.SupportsPositionalIndelScoresFor("B") {
		t.Errorf("Expected no scores for 'B'")
	}
	if exampleProfile.SupportsPositionalIndelScoresFor("A") {
		t.Errorf("Expected no support from this profile")
	}
}
