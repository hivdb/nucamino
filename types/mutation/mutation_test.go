package mutation

import (
	a "github.com/hivdb/nucamino/types/amino"
	c "github.com/hivdb/nucamino/types/codon"
	n "github.com/hivdb/nucamino/types/nucleic"
	"reflect"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestMakeMutation(t *testing.T) {
	result := MakeMutation(155, 797, []n.NucleicAcid{n.A, n.C, n.T}, a.S)
	expect := &Mutation{
		155, 797, "ACT", "T", &c.Codon{n.A, n.C, n.T}, "S", a.S,
		false, false, false, "...", "", "", nil,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{n.A, n.C, n.T}, a.T)
	if result != nil {
		t.Errorf(MSG_NOT_EQUAL, nil, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{n.A, n.T}, a.S)
	expect = &Mutation{
		155, 797, "A T", "NTSI", &c.Codon{n.A, n.N, n.T}, "S", a.S,
		false, false, true, ".-.", "", "", nil,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{}, a.S)
	expect = &Mutation{
		155, 797, "", "", nil, "S", a.S,
		false, true, false, "---", "", "", nil,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{n.A, n.C, n.T, n.A, n.C, n.T, n.R, n.C, n.T}, a.T)
	expect = &Mutation{
		155, 797, "ACT", "T", &c.Codon{n.A, n.C, n.T}, "T", a.T,
		true, false, false, ":::++++++", "ACTRCT", "T[TA]", []c.Codon{c.Codon{n.A, n.C, n.T}, c.Codon{n.R, n.C, n.T}},
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{n.T, n.A, n.R}, a.L)
	expect = &Mutation{
		155, 797, "TAR", "*", &c.Codon{n.T, n.A, n.R}, "L", a.L,
		false, false, false, "...", "", "", nil,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestGetInsertedCodons(t *testing.T) {
	result := MakeMutation(155, 797, []n.NucleicAcid{n.T, n.A, n.R}, a.L).GetInsertedCodons()
	if result != nil {
		t.Errorf(MSG_NOT_EQUAL, nil, result)
	}
}

func TestToString(t *testing.T) {
	result := MakeMutation(155, 797, []n.NucleicAcid{n.T, n.A, n.B}, a.L).ToString()
	expect := "L155Y*:TAB"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{n.T, n.B}, a.L).ToString()
	expect = "L155X:TB "
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{}, a.L).ToString()
	expect = "L155-"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeMutation(155, 797, []n.NucleicAcid{n.A, n.C, n.T, n.A, n.C, n.T, n.R, n.C, n.T}, a.T).ToString()
	expect = "T155T_T[TA]:ACT_ACTRCT"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}
