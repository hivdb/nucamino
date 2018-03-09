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
