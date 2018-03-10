package builtin

import (
	ap "github.com/hivdb/nucamino/alignmentprofile"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv1a"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hiv1b"
	"sort"
)

var profiles = map[string]ap.AlignmentProfile {
	"hcv1a": hcv1a.Profile,
	"hiv1b": hiv1b.Profile,
}

func Get(name string) (*ap.AlignmentProfile, bool) {
	profile, found := profiles[name]
	return &profile, found
}

func List() []string {
	keys := make([]string, 0, len(profiles))
	for k := range(profiles) {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}
