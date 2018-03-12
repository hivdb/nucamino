package cli

import "testing"

func TestValidOutputFormat(t *testing.T) {
	okCases := []string{"json", "tsv"}
	for _, c := range okCases {
		if !validOutputFormat(c) {
			t.Errorf("Expected %v to be a valid output format", c)
		}
	}
	errCases := []string{"csv", "xslx", "jpeg"}
	for _, c := range errCases {
		if validOutputFormat(c) {
			t.Errorf("Expected %v to not be a valid output format", c)
		}
	}
}
