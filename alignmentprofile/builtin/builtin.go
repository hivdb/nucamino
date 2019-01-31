package builtin

import (
	ap "github.com/hivdb/nucamino/alignmentprofile"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv1a"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv1b"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv2"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv3"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv4"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv5"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hcv6"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hiv1b"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hiv2a"
	"github.com/hivdb/nucamino/alignmentprofile/builtin/hiv2b"
	"sort"
)

var profiles = map[string]ap.AlignmentProfile{
	"hcv1a": hcv1a.Profile,
	"hcv1b": hcv1b.Profile,
	"hcv2":  hcv2.Profile,
	"hcv3":  hcv3.Profile,
	"hcv4":  hcv4.Profile,
	"hcv5":  hcv5.Profile,
	"hcv6":  hcv6.Profile,
	"hiv1b": hiv1b.Profile,
	"hiv2a": hiv2a.Profile,
	"hiv2b": hiv2b.Profile,
}

func Get(name string) (*ap.AlignmentProfile, bool) {
	profile, found := profiles[name]
	return &profile, found
}

func List() []string {
	keys := make([]string, 0, len(profiles))
	for k := range profiles {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}
