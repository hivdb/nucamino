package cmd

import (
	ap "github.com/hivdb/nucamino/alignmentprofile"
	"testing"
)

func TestGeneInGenes(t *testing.T) {
	var cases = []struct {
		gene     string
		genes    []ap.Gene
		expected bool
	}{
		{"a", []ap.Gene{"C", "B", "A"}, true},
		{"a", []ap.Gene{"C", "B", "A"}, true},
		{" a", []ap.Gene{"A", "B", "C"}, true},
		{"z", []ap.Gene{"A", "B", "C"}, false},
		{"z ", []ap.Gene{"A", "B", "C"}, false},
	}
	for _, cs := range cases {
		if geneInGenes(cs.gene, cs.genes) != cs.expected {
			t.Errorf("Expected %v in %v to be %v", cs.gene, cs.genes, cs.expected)
		}
	}
}

func TestAlignGetParameters(t *testing.T) {
	okCases := [][]string{
		[]string{"hcv1a", "ns3"},
		[]string{"hcv1a", "ns3, ns5b"},
		[]string{"hiv1b", "gag"},
		[]string{"hiv1b", "GAG,POL"},
	}
	for _, c := range okCases {
		_, _, err := alignGetParameters(c)
		if err != nil {
			t.Errorf("Expected alignGetParameters(...%v) to be okay", c)
		}
	}
	errCases := [][]string{
		// cobra ensures we get two args, so we only check cases with two args
		[]string{"asdf", "gag"},
		[]string{"hiv1b", "ns3"},
		[]string{"hcv1a", "n5a ns5b"},
	}
	for _, c := range errCases {
		_, _, err := alignGetParameters(c)
		if err == nil {
			t.Errorf("Expected alignGetParameters(...%v) to be an error", c)
		}
	}
}
