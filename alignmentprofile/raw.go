package alignmentprofile

// This file contains an intermediate representation of
// AlignmentProfile that's used during serialization and
// de-serialization. The serialization format is easier to read and
// write by hand; the AlignmentProfile structure is more suitable for
// calculation; the rawAlignmentProfile is both easy to
// serialize/deserialize and easy to convert into an AlignmentProfile.

import (
	"fmt"
	a "github.com/hivdb/nucamino/types/amino"
)

// This structure is a de-serialization target that the YAML package
// uses to parse the GeneIndelScores in a serialized profile.
type rawIndelScore struct {
	Kind     string
	Position int
	Open     int
	Extend   int
}

func (t *rawIndelScore) UnmarshalYAML(unmarshal func(interface{}) error) error {
	bucket := make([]interface{}, 4)
	e := unmarshal(&bucket)
	t.Kind = bucket[0].(string)
	t.Position = bucket[1].(int)
	t.Open = bucket[2].(int)
	t.Extend = bucket[3].(int)
	if e != nil {
		return e
	}
	return nil
}

// This type alias lets us implement the sorting interface for
// []rawIndelScore. We sort it first by position (with lower positions
// first) and then by kind (with insertions before deletions).
type byPositionAndKind []rawIndelScore

func (a byPositionAndKind) Len() int {
	return len(a)
}

func (a byPositionAndKind) Swap(i, j int) {
	a[i], a[j] = a[j], a[i]
}

func (a byPositionAndKind) Less(i, j int) bool {
	if a[i].Position < a[j].Position {
		return true
	} else if a[i].Position == a[j].Position {
		return a[i].Kind == "ins" && a[j].Kind == "del"
	}
	return false
}

// This is an intermediate datatype between an AlignmentProfile and
// the YAML that represents it. The YAML is formatted for editing,
// while the AlignmentProfile is formatted for ease of
// calculation. This structure can be deserialized form YAML,
// converted to an AlignmentProfile, or contructed from an
// AlignmentProfile.
type rawAlignmentProfile struct {
	StopCodonPenalty         int                        `yaml:"StopCodonPenalty"`
	GapOpeningPenalty        int                        `yaml:"GapOpeningPenalty"`
	GapExtensionPenalty      int                        `yaml:"GapExtensionPenalty"`
	IndelCodonOpeningBonus   int                        `yaml:"IndelCodonOpeningBonus"`
	IndelCodonExtensionBonus int                        `yaml:"IndelCodonExtensionBonus"`
	RawIndelScores           map[string][]rawIndelScore `yaml:"PositionalIndelScores,flow"`
	ReferenceSequences       map[string]string          `yaml:"ReferenceSequences"`
}

// Construct a GenePositionalIndelScores instance from a
// rawAlignentProfile (which has presumably been parsed from YAML).
func (rawProfile rawAlignmentProfile) geneIndelScores() (*GenePositionalIndelScores, error) {
	geneIndelScores := make(GenePositionalIndelScores)
	for geneSrc, rawIndelScores := range rawProfile.RawIndelScores {
		indelScores := make(PositionalIndelScores)
		for _, indelScore := range rawIndelScores {
			var scoreKeySign int
			if indelScore.Kind == "ins" {
				scoreKeySign = 1
			} else if indelScore.Kind == "del" {
				scoreKeySign = -1
			} else {
				msgFmt := "Unknown indel score kind '%v' (expecting 'ins' or 'del')"
				return nil, fmt.Errorf(msgFmt, indelScore.Kind)
			}
			indelKey := scoreKeySign * indelScore.Position
			indelScores[indelKey] = [2]int{indelScore.Open, indelScore.Extend}
		}
		geneIndelScores[Gene(geneSrc)] = indelScores
	}
	return &geneIndelScores, nil
}

// Construct an AlignmentProfile from a rawAlignmentProfile
func (raw rawAlignmentProfile) asProfile() (*AlignmentProfile, error) {
	var profile AlignmentProfile
	profile.StopCodonPenalty = raw.StopCodonPenalty
	profile.GapOpeningPenalty = raw.GapOpeningPenalty
	profile.GapExtensionPenalty = raw.GapExtensionPenalty
	profile.IndelCodonOpeningBonus = raw.IndelCodonOpeningBonus
	profile.IndelCodonExtensionBonus = raw.IndelCodonExtensionBonus

	if len(raw.ReferenceSequences) == 0 {
		return nil, fmt.Errorf("Missing key: ReferenceSequences")
	} else {
		profile.ReferenceSequences = make(ReferenceSeqs)
		for geneSrc, aaSrc := range raw.ReferenceSequences {
			profile.ReferenceSequences[Gene(geneSrc)] = a.ReadString(aaSrc)
		}
	}

	if raw.RawIndelScores == nil || len(raw.RawIndelScores) == 0 {
		profile.GeneIndelScores = nil
	} else {
		geneIndelScores, err := raw.geneIndelScores()
		if err != nil {
			return nil, err
		}
		profile.GeneIndelScores = *geneIndelScores
	}

	return &profile, nil
}
