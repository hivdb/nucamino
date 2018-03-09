package alignmentprofile

import (
	"fmt"
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
	RawIndelScores           map[string][]rawIndelScore `yaml:"GeneIndelScores,flow"`
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
