package alignment

import (
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	// c "github.com/hivdb/nucamino/types/codon"
	n "github.com/hivdb/nucamino/types/nucleic"
	"reflect"
	"testing"
)

const (
	MSG_NOT_EQUAL = "Expect %#v but received %#v"
)

var (
	NSEQ = n.ReadString("ACAGTRTTAGTAGGACCTACACCTGCCAACATAATTGGAAGAAATCTGTTGACYCAG")
	ASEQ = a.ReadString("TVLVGPTPVNIIGRNLLTQ")
)

func TestNewAlignment(t *testing.T) {
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	result := NewAlignment(NSEQ, ASEQ, handler)
	expect := &Alignment{
		q: -1000, r: -200,
		nSeq: NSEQ, aSeq: ASEQ, nSeqLen: 57, aSeqLen: 19,
		scoreHandler: handler, nwMatrix: []int{},
		endPosN: 57, endPosA: 19, maxScore: 9100,
		isSimpleAlignment:             true,
		supportPositionalIndel:        false,
		constIndelCodonOpeningScore:   0,
		constIndelCodonExtensionScore: 200,
	}
	result.report = nil
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestGetMatrixIndex(t *testing.T) {
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	aln := NewAlignment(NSEQ, ASEQ, handler)
	result := aln.getMatrixIndex(INS, 14, 5)
	expect := 1445
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = aln.getMatrixIndex(DEL, 14, 5)
	expect = 2605
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
	result = aln.getMatrixIndex(GENERAL, 14, 5)
	expect = 285
	if result != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestGetTypedPos(t *testing.T) {
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	aln := NewAlignment(NSEQ, ASEQ, handler)
	st, pn, pa := aln.getTypedPos(1445)
	if st != INS || pn != 14 || pa != 5 {
		t.Errorf(MSG_NOT_EQUAL, []int{int(INS), 14, 5}, []int{int(st), pn, pa})
	}
	st, pn, pa = aln.getTypedPos(2605)
	if st != DEL || pn != 14 || pa != 5 {
		t.Errorf(MSG_NOT_EQUAL, []int{int(DEL), 14, 5}, []int{int(st), pn, pa})
	}
	st, pn, pa = aln.getTypedPos(285)
	if st != GENERAL || pn != 14 || pa != 5 {
		t.Errorf(MSG_NOT_EQUAL, []int{int(GENERAL), 14, 5}, []int{int(st), pn, pa})
	}
}
