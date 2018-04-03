package alignmentprofile

import "strings"

func (g Gene) Matches(s string) bool {
	candidate := strings.ToUpper(strings.TrimSpace(s))
	return candidate == string(g)
}
