package utils

import "testing"

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestStripWhiteSpace(t *testing.T) {
	result := StripWhiteSpace("   AB  C\t\n\r\n  ")
	expectResult := "ABC"
	if result != expectResult {
		t.Errorf(MSG_NOT_EQUAL, expectResult, result)
	}

}
