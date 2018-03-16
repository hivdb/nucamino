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
	StopCodonPenalty:         0,
	GapOpeningPenalty:        1,
	GapExtensionPenalty:      2,
	IndelCodonOpeningBonus:   3,
	IndelCodonExtensionBonus: 4,
	RawIndelScores: map[string][]rawIndelScore{
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
	ReferenceSequences: map[string]string{
		"A": "SGSWLRD",
		"B": "CPPDSDVESYSSMPPLEGEPGDPD",
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
	RawIndelScores: map[string][]rawIndelScore{
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

func TestRawIndelScoreSortOrder(t *testing.T) {
	lessCases := [][]rawIndelScore{
		{
			{Position: 1},
			{Position: 2},
		},
		{
			{Position: 1, Kind: "ins"},
			{Position: 1, Kind: "del"},
		},
		{
			{Position: 1, Kind: "del"},
			{Position: 2, Kind: "ins"},
		},
	}
	for _, c := range lessCases {
		if !byPositionAndKind(c).Less(0, 1) {
			t.Errorf("Expected %v to be less than %v", c[0], c[1])
		}
	}
	notLessCases := [][]rawIndelScore{
		{
			{Position: 2},
			{Position: 1},
		},
		{
			{Position: 1, Kind: "del"},
			{Position: 1, Kind: "ins"},
		},
	}
	for _, c := range notLessCases {
		if byPositionAndKind(c).Less(0, 1) {
			t.Errorf("Expected %v to be greater than %v", c[0], c[1])
		}
	}
}

func TestRoundTipToAlignmentProfile(t *testing.T) {
	profile, err := exampleRawProfile.asProfile()
	if err != nil {
		t.Errorf("Unexpected error while converting raw profile to profile: %v", err)
		t.FailNow()
	}
	constructed := profile.asRaw()
	if !reflect.DeepEqual(constructed, exampleRawProfile) {
		t.Errorf("%v != %v", constructed, exampleRawProfile)
	}
}
