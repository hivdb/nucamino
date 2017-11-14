package alignment

import (
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	f "github.com/hivdb/nucamino/types/frameshift"
	m "github.com/hivdb/nucamino/types/mutation"
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

func TestNewAlignmentWithIns(t *testing.T) {
	NSEQ_INS := n.ReadString("ACAGTRTTAGTAGGACCTTTTACACCTGCCAACATAATTGGAAGAAATCTGTTGACYCAG")
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	result := NewAlignment(NSEQ_INS, ASEQ, handler)
	expect := &Alignment{
		q: -1000, r: -200,
		nSeq: NSEQ_INS, aSeq: ASEQ, nSeqLen: 60, aSeqLen: 19,
		scoreHandler: handler, nwMatrix: []int{},
		endPosN: 60, endPosA: 19, maxScore: 7700,
		isSimpleAlignment:             false,
		supportPositionalIndel:        false,
		constIndelCodonOpeningScore:   0,
		constIndelCodonExtensionScore: 200,
	}
	result.nwMatrix = []int{}
	result.report = nil
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestNewAlignmentWithInsFs(t *testing.T) {
	NSEQ_INSFS := n.ReadString("ACAGTRTTAGTAGGACCTTTACACCTGCCAACATAATTGGAAGAAATCTGTTGACYCAG")
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	result := NewAlignment(NSEQ_INSFS, ASEQ, handler)
	expect := &Alignment{
		q: -1000, r: -200,
		nSeq: NSEQ_INSFS, aSeq: ASEQ, nSeqLen: 59, aSeqLen: 19,
		scoreHandler: handler, nwMatrix: []int{},
		endPosN: 59, endPosA: 19, maxScore: 7700,
		isSimpleAlignment:             false,
		supportPositionalIndel:        false,
		constIndelCodonOpeningScore:   0,
		constIndelCodonExtensionScore: 200,
	}
	result.nwMatrix = []int{}
	result.report = nil
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestNewAlignmentWithDelFs(t *testing.T) {
	NSEQ_DELFS := n.ReadString("AAGTRTTAGTAGGACCTACACCTGCCAACATAATTGGAGAAATCTGTTGACYCAG")
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	result := NewAlignment(NSEQ_DELFS, ASEQ, handler)
	expect := &Alignment{
		q: -1000, r: -200,
		nSeq: NSEQ_DELFS, aSeq: ASEQ, nSeqLen: 55, aSeqLen: 19,
		scoreHandler: handler, nwMatrix: []int{},
		endPosN: 55, endPosA: 19, maxScore: 7000,
		isSimpleAlignment:             false,
		supportPositionalIndel:        false,
		constIndelCodonOpeningScore:   0,
		constIndelCodonExtensionScore: 200,
	}
	result.nwMatrix = []int{}
	result.report = nil
	if !reflect.DeepEqual(expect, result) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}

func TestNewAlignmentWithDel(t *testing.T) {
	NSEQ_DEL := n.ReadString("ACAGTRTTAGTAGGACCTACACCTAACATAATTGGAAGAAATCTGTTGACYCA")
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	result := NewAlignment(NSEQ_DEL, ASEQ, handler)
	expect := &Alignment{
		q: -1000, r: -200,
		nSeq: NSEQ_DEL, aSeq: ASEQ, nSeqLen: 53, aSeqLen: 19,
		scoreHandler: handler, nwMatrix: []int{},
		endPosN: 53, endPosA: 19, maxScore: 7450,
		isSimpleAlignment:             false,
		supportPositionalIndel:        false,
		constIndelCodonOpeningScore:   0,
		constIndelCodonExtensionScore: 200,
	}
	result.nwMatrix = []int{}
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

func TestGetReport(t *testing.T) {
	nseq := n.ReadString("ACAGTRTTAGTAGGACCTACACCTAACATAATTGGAAGAAAAAATCTGTTGACYCA")
	handler := h.New(4, 10, 2, 0, 2, nil, false)
	aln := NewAlignment(nseq, ASEQ, handler)
	result := aln.GetReport()
	expect := &AlignmentReport{
		FirstAA: 1,
		FirstNA: 1,
		LastAA:  18,
		LastNA:  54,
		Mutations: []m.Mutation{
			*m.MakeMutation(9, 25, []n.NucleicAcid{}, a.V),
			*m.MakeMutation(14, 37, []n.NucleicAcid{n.A, n.G, n.A, n.A, n.A, n.A}, a.R),
		},
		FrameShifts: []f.FrameShift{},
		AlignedSites: []AlignedSite{
			AlignedSite{1, 1, 3},
			AlignedSite{2, 4, 3},
			AlignedSite{3, 7, 3},
			AlignedSite{4, 10, 3},
			AlignedSite{5, 13, 3},
			AlignedSite{6, 16, 3},
			AlignedSite{7, 19, 3},
			AlignedSite{8, 22, 3},
			AlignedSite{9, 25, 0},
			AlignedSite{10, 25, 3},
			AlignedSite{11, 28, 3},
			AlignedSite{12, 31, 3},
			AlignedSite{13, 34, 3},
			AlignedSite{14, 37, 6},
			AlignedSite{15, 43, 3},
			AlignedSite{16, 46, 3},
			AlignedSite{17, 49, 3},
			AlignedSite{18, 52, 3},
		},
		AminoAcidsLine:    "T  V  L  V  G  P  T  P  V  N  I  I  G  R     N  L  L  T  ",
		ControlLine:       "::::::::::::::::::::::::---:::::::::::::::+++::::::::::::",
		NucleicAcidsLine:  "ACAGTRTTAGTAGGACCTACACCT   AACATAATTGGAAGAAAAAATCTGTTGACY",
		IsSimpleAlignment: false,
	}
	if !reflect.DeepEqual(result, expect) {
		t.Errorf(MSG_NOT_EQUAL, expect, result)
	}
}
