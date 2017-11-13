package nucleic

import (
	"reflect"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestToString(t *testing.T) {
	if A.ToString() != "A" {
		t.Errorf(MSG_NOT_EQUAL, "A", A.ToString())
	}
	if B.ToString() != "B" {
		t.Errorf(MSG_NOT_EQUAL, "B", B.ToString())
	}
}

func TestIsAmbiguous(t *testing.T) {
	if A.IsAmbiguous() {
		t.Errorf(MSG_NOT_EQUAL, false, true)
	}
	if T.IsAmbiguous() {
		t.Errorf(MSG_NOT_EQUAL, false, true)
	}
	if !W.IsAmbiguous() {
		t.Errorf(MSG_NOT_EQUAL, true, false)
	}
	if !N.IsAmbiguous() {
		t.Errorf(MSG_NOT_EQUAL, true, false)
	}
}

func TestReadString(t *testing.T) {
	result := ReadString("ACGTwWSM ikKRYBDHVN")
	expect := []NucleicAcid{A, C, G, T, W, W, S, M, N, K, K, R, Y, B, D, H, V, N}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestWriteString(t *testing.T) {
	result := WriteString([]NucleicAcid{A, C, G, T, W, W, S, M, N, K, K, R, Y, B, D, H, V, N})
	expect := "ACGTWWSMNKKRYBDHVN"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestGetUnambiguousNucleicAcids(t *testing.T) {
	result := GetUnambiguousNucleicAcids(A)
	expect := []NucleicAcid{A}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = GetUnambiguousNucleicAcids(R)
	expect = []NucleicAcid{A, G}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = GetUnambiguousNucleicAcids(B)
	expect = []NucleicAcid{C, G, T}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = GetUnambiguousNucleicAcids(N)
	expect = []NucleicAcid{A, C, G, T}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}
