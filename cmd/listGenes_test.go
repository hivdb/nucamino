package cmd

import "testing"

func TestListGenesErrorHandling(t *testing.T) {
	var builtinProfile = "hcv1a"
	var err error

	err = listGenes(nil, []string{"not-a-profile"})
	if err == nil {
		t.Errorf("Expected missing gene-name to raise an error")
	}

	var invalidRegex = "a[bc"
	err = listGenes(nil, []string{builtinProfile, invalidRegex})
	if err == nil {
		t.Errorf("Expected '%v' to be an invalid regex", invalidRegex)
	}

	err = listGenes(nil, []string{builtinProfile})
	if err != nil {
		t.Fail()
	}

	var validRegex = "3$"
	err = listGenes(nil, []string{builtinProfile, validRegex})
	if err != nil {
		t.Fail()
	}
}
