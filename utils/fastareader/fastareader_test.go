package fastareader_test

import (
	"fmt"
	n "github.com/hivdb/nucamino/types/nucleic"
	"github.com/hivdb/nucamino/utils/fastareader"
	"strings"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

func TestIterSequences(t *testing.T) {
	reader := strings.NewReader(`
>TestSeq1
ACGT
TG CA
;comment1
>TestSeq2
;comment2
WSMK.
  BDH-VN`)
	seqChan := fastareader.IterSequences(reader)
	seq := <-seqChan
	expectName := "TestSeq1"
	if seq.Name != expectName {
		t.Errorf(MSG_NOT_EQUAL, expectName, seq.Name)
	}
	expectSeq := fmt.Sprintf(
		"%#v", []n.NucleicAcid{n.A, n.C, n.G, n.T, n.T, n.G, n.C, n.A})
	if fmt.Sprintf("%#v", seq.Sequence) != expectSeq {
		t.Errorf(MSG_NOT_EQUAL, expectSeq, seq.Sequence)
	}
	seq = <-seqChan
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
