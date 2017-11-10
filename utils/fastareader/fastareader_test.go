package fastareader

import (
	"fmt"
	n "github.com/hivdb/nucamino/types/nucleic"
	"strings"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestReadSequences(t *testing.T) {
	reader := strings.NewReader(`
>TestSeq1
ACGT
TG CA
;comment1
>TestSeq2
;comment2
WSMK.
  BDH-VN`)
	seqs := ReadSequences(reader)
	seq := seqs[0]
	expectName := "TestSeq1"
	if seq.Name != expectName {
		t.Errorf(MSG_NOT_EQUAL, expectName, seq.Name)
	}
	expectSeq := fmt.Sprintf(
		"%#v", []n.NucleicAcid{n.A, n.C, n.G, n.T, n.T, n.G, n.C, n.A})
	if fmt.Sprintf("%#v", seq.Sequence) != expectSeq {
		t.Errorf(MSG_NOT_EQUAL, expectSeq, seq.Sequence)
	}
	seq = seqs[1]
	expectName = "TestSeq2"
	if seq.Name != expectName {
		t.Errorf(MSG_NOT_EQUAL, expectName, seq.Name)
	}
	expectSeq = fmt.Sprintf(
		"%#v", []n.NucleicAcid{
			n.W, n.S, n.M, n.K, n.N /*.*/, n.B, n.D, n.H, n.N /*-*/, n.V, n.N})
	if fmt.Sprintf("%#v", seq.Sequence) != expectSeq {
		t.Errorf(MSG_NOT_EQUAL, expectSeq, seq.Sequence)
	}
}

func TestReadSequencesBranch1(t *testing.T) {
	reader := strings.NewReader(`
>
ACGT
TG CA
>
AAAA`)
	seqs := ReadSequences(reader)
	seq := seqs[0]
	expectName := "unnamed sequence 1"
	if seq.Name != expectName {
		t.Errorf(MSG_NOT_EQUAL, expectName, seq.Name)
	}
	seq = seqs[1]
	expectName = "unnamed sequence 2"
	if seq.Name != expectName {
		t.Errorf(MSG_NOT_EQUAL, expectName, seq.Name)
	}
}

func TestReadSequencesBranch2(t *testing.T) {
	reader := strings.NewReader(`
ACGT
TG CA
AAAA`)
	seqs := ReadSequences(reader)
	seq := seqs[0]
	expectName := "unnamed sequence"
	if seq.Name != expectName {
		t.Errorf(MSG_NOT_EQUAL, expectName, seq.Name)
	}
}
