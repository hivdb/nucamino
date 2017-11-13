package codon

import (
	a "github.com/hivdb/nucamino/types/amino"
	. "github.com/hivdb/nucamino/types/nucleic"
	"reflect"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestHasCommonCodon(t *testing.T) {
	if !hasCommonCodon(
		[]Codon{
			Codon{T, G, T},
			Codon{A, G, T},
		},
		[]Codon{
			Codon{G, G, T},
			Codon{G, G, C},
			Codon{A, G, T},
		},
	) {
		t.Errorf(MSG_NOT_EQUAL, true, false)
	}
	if hasCommonCodon(
		[]Codon{},
		[]Codon{
			Codon{G, G, T},
		},
	) {
		t.Errorf(MSG_NOT_EQUAL, false, true)
	}
	if hasCommonCodon(
		[]Codon{
			Codon{A, G, T},
		},
		[]Codon{},
	) {
		t.Errorf(MSG_NOT_EQUAL, false, true)
	}
	if hasCommonCodon(
		[]Codon{
			Codon{C, T, T},
			Codon{C, T, G},
			Codon{A, G, T},
		},
		[]Codon{
			Codon{T, C, T},
			Codon{G, G, T},
		},
	) {
		t.Errorf(MSG_NOT_EQUAL, false, true)
	}
}

func TestFindBestMatch(t *testing.T) {
	result := FindBestMatch([]NucleicAcid{C}, a.T)
	expect := Codon{N, C, N}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{G}, a.Y)
	expect = Codon{G, N, N}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{C}, a.L)
	expect = Codon{C, N, N}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{B}, a.I)
	expect = Codon{N, B, N}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{G, C}, a.Y)
	expect = Codon{G, N, C}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{C, A}, a.T)
	expect = Codon{N, C, A}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{G, A, A}, a.L)
	expect = Codon{G, A, A}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{T, T, A}, a.F)
	expect = Codon{T, T, A}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{C, R}, a.P)
	expect = Codon{N, C, R}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{C, R}, a.S)
	expect = Codon{N, C, R}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = FindBestMatch([]NucleicAcid{R, V}, a.V)
	expect = Codon{R, N, V}
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestIsStopCodon(t *testing.T) {
	if (&Codon{T, T, T}).IsStopCodon() {
		t.Errorf(MSG_NOT_EQUAL, false, true)
	}
	if !(&Codon{T, A, A}).IsStopCodon() {
		t.Errorf(MSG_NOT_EQUAL, true, false)
	}
	if !(&Codon{T, A, R}).IsStopCodon() {
		t.Errorf(MSG_NOT_EQUAL, true, false)
	}
}

func TestToAminoAcidsText(t *testing.T) {
	result := (&Codon{A, A, A}).ToAminoAcidsText()
	expect := "K"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = (&Codon{T, A, A}).ToAminoAcidsText()
	expect = "*"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = (&Codon{R, T, T}).ToAminoAcidsText()
	expect = "IV"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = (&Codon{R, T, R}).ToAminoAcidsText()
	expect = "IMV"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestIsAmbiguous(t *testing.T) {
	if (&Codon{A, A, A}).IsAmbiguous() {
		t.Errorf(MSG_NOT_EQUAL, false, true)
	}
	if !(&Codon{A, A, R}).IsAmbiguous() {
		t.Errorf(MSG_NOT_EQUAL, true, false)
	}
	if !(&Codon{N, A, A}).IsAmbiguous() {
		t.Errorf(MSG_NOT_EQUAL, true, false)
	}
}

func TestToString(t *testing.T) {
	result := (&Codon{A, A, A}).ToString()
	expect := "AAA"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = (&Codon{R, A, N}).ToString()
	expect = "RAN"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}
