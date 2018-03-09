package alignmentprofile

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
	GeneIndelScores          map[string][]rawIndelScore `yaml:"GeneIndelScore,flow"`
	ReferenceSequences       map[string]string          `yaml:"ReferenceSequences"`
}
