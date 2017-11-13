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

func TestPadRightSpace(t *testing.T) {
	result := PadRightSpace("AB ", 8)
	expect := "AB      "
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}
