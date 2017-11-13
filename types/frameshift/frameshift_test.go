package frameshift

import (
	n "github.com/hivdb/nucamino/types/nucleic"
	"reflect"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestNew(t *testing.T) {
	result := New(155, 677, []n.NucleicAcid{}, DELETION, 2)
	expect := &FrameShift{
		155,
		677,
		[]n.NucleicAcid{},
		"",
		DELETION,
		false,
		true,
		2,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = New(155, 677, []n.NucleicAcid{n.A, n.R}, INSERTION, 2)
	expect = &FrameShift{
		155,
		677,
		[]n.NucleicAcid{n.A, n.R},
		"AR",
		INSERTION,
		true,
		false,
		2,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestMakeFrameShift(t *testing.T) {
	result := MakeFrameShift(155, 677, []n.NucleicAcid{n.A, n.R, n.G})
	if result != nil {
		t.Errorf(MSG_NOT_EQUAL, nil, result)
	}
	result = MakeFrameShift(155, 677, []n.NucleicAcid{n.A, n.R, n.G, n.T})
	expect := &FrameShift{
		155,
		677 + 4/3*3,
		[]n.NucleicAcid{n.T},
		"T",
		INSERTION,
		true,
		false,
		1,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeFrameShift(155, 677, []n.NucleicAcid{n.A, n.R})
	expect = &FrameShift{
		155,
		677 + 2,
		[]n.NucleicAcid{},
		"",
		DELETION,
		false,
		true,
		1,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestToString(t *testing.T) {
	result := MakeFrameShift(155, 677, []n.NucleicAcid{n.A, n.R, n.G, n.T}).ToString()
	expect := "155ins1bp_T"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = MakeFrameShift(155, 677, []n.NucleicAcid{n.A, n.R}).ToString()
	expect = "155del1bp"
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}
