package alignmentprofile

import "testing"

func TestUnmarshalEmpty(t *testing.T) {
	src := ""
	_, err := Parse(src)
	if err == nil {
		t.Errorf("Expected error on empty YAML file")
	}
}

func TestInvalidProfile(t *testing.T) {
	src := `StopCodonPenalty: 0
GapOpeningPenalty: 0
GapExtensionPenalty: 0
IndelCodonOpeningBonus: 0
IndelCodonExtensionBonus: 0`
	_, err := Parse(src)
	if err == nil {
		t.Errorf("Expected error when missing ReferenceSequences")
	}
}
