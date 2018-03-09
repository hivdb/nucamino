package alignmentprofile

import (
	yaml "gopkg.in/yaml.v2"
	"reflect"
	"testing"
)

func TestUnmarshalRawIndelScore(t *testing.T) {
	src := []byte("[ins, 0, 1, 2]")
	expected := rawIndelScore{
		Kind:     "ins",
		Position: 0,
		Open:     1,
		Extend:   2,
	}
	var constructed rawIndelScore
	err := yaml.Unmarshal(src, &constructed)
	if err != nil {
		t.Errorf("Unexpected error: %v", err)
	}
	if !reflect.DeepEqual(constructed, expected) {
		t.Errorf("%+v != %+v", constructed, expected)
	}
}

var exampleRawProfile = rawAlignmentProfile{
	GeneIndelScores: map[string][]rawIndelScore{
		"A": []rawIndelScore{
			rawIndelScore{
				Kind:     "ins",
				Position: 1,
				Open:     2,
				Extend:   3,
			},
			rawIndelScore{
				Kind:     "del",
				Position: 4,
				Open:     5,
				Extend:   6,
			},
		},
	},
}

func TestGeneIndelScoresFromRawProfile(t *testing.T) {
	constructed, err := exampleRawProfile.geneIndelScores()
	expected := &GenePositionalIndelScores{
		Gene("A"): PositionalIndelScores{
			1:  [2]int{2, 3},
			-4: [2]int{5, 6},
		},
	}
	if err != nil {
		t.Errorf("Unexpected Error: %v", err)
	}
	if !reflect.DeepEqual(constructed, expected) {
		t.Errorf("%+v != %+v", constructed, expected)
	}
}

var exampleInvalidRawProfile = rawAlignmentProfile{
	GeneIndelScores: map[string][]rawIndelScore{
		"A": []rawIndelScore{
			rawIndelScore{
				Kind:     "foo",
				Position: 1,
				Open:     2,
				Extend:   3,
			},
		},
	},
}

func TestGeneIndelScoresFromRawProfileWithInvalidScoreKind(t *testing.T) {
	_, err := exampleInvalidRawProfile.geneIndelScores()
	if err == nil {
		t.Errorf("Expected an error on indel kind 'foo'")
	}
}
