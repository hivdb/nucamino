package alignmentprofile

import (
	a "github.com/hivdb/nucamino/types/amino"
	"testing"
)

var exampleProfile = AlignmentProfile{
	ReferenceSequences: ReferenceSeqs{
		"A": a.ReadString("T"),
		"B": a.ReadString("V"),
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
