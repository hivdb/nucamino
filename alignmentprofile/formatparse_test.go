package alignmentprofile

import (
	"testing"
	"reflect"
)


func TestFormatParseRoundTrip(t *testing.T) {
	formatted := Format(exampleProfile)
	parsed, err := Parse(formatted)
	if err != nil {
		t.Errorf("Unexpected error while parsing formatted profile: %v", err)
		t.FailNow()
	}
	if !reflect.DeepEqual(*parsed, exampleProfile) {
		t.Errorf("%v != %v", *parsed, exampleProfile)
	}
}

var exampleProfileYAML =`StopCodonPenalty: 1
GapOpeningPenalty: 2
GapExtensionPenalty: 3
IndelCodonOpeningBonus: 4
IndelCodonExtensionBonus: 5
ReferenceSequences:
  A:
    TTALIEPPVYPIVEHSDEKTAHEEH
  B:
    CSNELVISHEADPVWRSAVLRGAP
PositionalIndelScores:
  A:
    - [ ins, 3, 4, 5 ]
    - [ ins, 6, 7, 8 ]
    - [ del, 6, 7, 8 ]
    - [ del, 9, 10, 11 ]
  B:
    - [ ins, 2, 1, 2 ]
`

func TestParseFormatRoundTrip(t *testing.T) {
	parsed, err := Parse(exampleProfileYAML)
	if err != nil {
		t.Errorf("Unexpected error while parsing example YAML: %v", err)
		t.FailNow()
	}
	formatted := Format(*parsed)
	if formatted != exampleProfileYAML {
		t.Errorf("%v != %v", formatted, exampleProfileYAML)
	}
}
