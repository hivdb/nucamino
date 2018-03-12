package alignmentprofile

import "testing"

func TestGeneMatching(t *testing.T) {
	matchingCases := []struct {
		g Gene
		s string
	}{
		{Gene("A"), "a"},
		{Gene("A"), "A"},
		{Gene("A"), " a "},
		{Gene("A"), "A "},
	}
	for _, c := range matchingCases {
		if !c.g.Matches(c.s) {
			t.Errorf("%v doesn't match %v", c.g, c.s)
		}
	}
	mismatchCases := []struct {
		g Gene
		s string
	}{
		{Gene("A"), "B"},
		{Gene("A"), "b"},
		{Gene("A"), "A1"},
		{Gene("A"), " A1"},
		{Gene("A"), " A 1"},
	}
	for _, c := range mismatchCases {
		if c.g.Matches(c.s) {
			t.Errorf("%v matches %v", c.g, c.s)
		}
	}
}
