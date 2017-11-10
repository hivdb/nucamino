package amino

import (
	"reflect"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestToString(t *testing.T) {
	aminoAcids := []AminoAcid{
		A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
	}
	expects := []string{
		"A", "C", "D", "E", "F", "G", "H", "I", "K",
		"L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
	}
	for i := 0; i < len(aminoAcids); i++ {
		aaText := ToString(aminoAcids[i])
		if aaText != expects[i] {
			t.Errorf("Expect string %#v but got %#v", expects[i], aaText)
		}
	}
}

func TestReadString(t *testing.T) {
	result := ReadString("ABCDEFG")
	expect := []AminoAcid{A, C, D, E, F, G} // ignore "B"
	if !reflect.DeepEqual(result, expect) {
		t.Errorf("Expect %#v but got %#v", expect, result)
	}
}

func TestWriteString(t *testing.T) {
	result := WriteString([]AminoAcid{A, C, D, E, F, G, V, W, Y, Y})
	expect := "ACDEFGVWYY"
	if result != expect {
		t.Errorf("Expect %#v but got %#v", expect, result)
	}
}
