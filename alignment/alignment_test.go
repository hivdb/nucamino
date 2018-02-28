package alignment

import (
	"errors"
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
	NSEQ                 = n.ReadString("ACAGTRTTAGTAGGACCTACACCTGCCAACATAATTGGAAGAAATCTGTTGACYCAG")
	ASEQ                 = a.ReadString("TVLVGPTPVNIIGRNLLTQ")
	BASIC_HANDLER_PARAMS = h.GeneralScoreHandlerParams{
		StopCodonPenalty:              4,
		GapOpenPenalty:                10,
		GapExtensionPenalty:           2,
		IndelCodonOpeningBonus:        0,
		IndelCodonExtensionBonus:      2,
		PositionalIndelScores:         nil,
		SupportsPositionalIndelScores: false,
	}
)

func TestNewAlignment(t *testing.T) {
	handler := h.New(BASIC_HANDLER_PARAMS)
	result, _ := NewAlignment(NSEQ, ASEQ, handler)
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
	handler := h.New(BASIC_HANDLER_PARAMS)
	result, _ := NewAlignment(NSEQ_INS, ASEQ, handler)
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
	handler := h.New(BASIC_HANDLER_PARAMS)
	result, _ := NewAlignment(NSEQ_INSFS, ASEQ, handler)
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
	handler := h.New(BASIC_HANDLER_PARAMS)
	result, _ := NewAlignment(NSEQ_DELFS, ASEQ, handler)
	expect := &Alignment{
		q: -1000, r: -200,
		nSeq: NSEQ_DELFS, aSeq: ASEQ, nSeqLen: 55, aSeqLen: 19,
		scoreHandler: handler, nwMatrix: []int{},
		endPosN: 55, endPosA: 19, maxScore: 6500,
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
	NSEQ_DEL := n.ReadString("ACAGTRTTAGTAGGACCTACACCT   AACATAATTGGAAGAAATCTGTTGACYCA")
	handler := h.New(BASIC_HANDLER_PARAMS)
	result, _ := NewAlignment(NSEQ_DEL, ASEQ, handler)
	expect := &Alignment{
		q: -1000, r: -200,
		nSeq: NSEQ_DEL, aSeq: ASEQ, nSeqLen: 53, aSeqLen: 19,
		scoreHandler: handler, nwMatrix: []int{},
		endPosN: 53, endPosA: 19, maxScore: 7200,
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
	handler := h.New(BASIC_HANDLER_PARAMS)
	aln, _ := NewAlignment(NSEQ, ASEQ, handler)
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
	handler := h.New(BASIC_HANDLER_PARAMS)
	aln, _ := NewAlignment(NSEQ, ASEQ, handler)
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

func TestMisaligned(t *testing.T) {
	nseq := n.ReadString("AAAAAAAAAAAAAAAAAA")
	aseq := a.ReadString("GGGGGGGGGGGGGGGGGG")
	handler := h.New(BASIC_HANDLER_PARAMS)
	aln, err := NewAlignment(nseq, aseq, handler)
	if aln != nil {
		t.Errorf(MSG_NOT_EQUAL, nil, aln)
	}
	errExpect := errors.New("sequence misaligned")
	if !reflect.DeepEqual(err, errExpect) {
		t.Errorf(MSG_NOT_EQUAL, errExpect, err)
	}
}

func TestGetReport(t *testing.T) {
	nseq := n.ReadString("ACAGTRTTAGTAGGACCTACACCTAACATAATTGGAAGAAAAAATCTGTTGACYCA")
	handler := h.New(BASIC_HANDLER_PARAMS)
	aln, _ := NewAlignment(nseq, ASEQ, handler)
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

func TestPosInsCodonScore(t *testing.T) {
	nseq := n.ReadString("ACAGTRTTAGTAGGACCTACACCTttttttGCCAACATAATTGGAAGAAATCTGTTGACYCAG")
	handler := h.New(BASIC_HANDLER_PARAMS)
	aln, _ := NewAlignment(nseq, ASEQ, handler)
	result := aln.GetReport()
	//           1  2  3  4  5  6  7  8        9 10 11 12 13 14 15 16 17 18 19
	//         ACAGTRTTAGTAGGACCTACACCTttttttGCCAACATAATTGGAAGAAATCTGTTGACYCAG
	expect := "::::::::::::::::::::::::++++++...::::::::::::::::::::::::::::::"
	if result.ControlLine != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result.ControlLine)
	}
	handler = h.New(h.GeneralScoreHandlerParams{
		StopCodonPenalty:         4,
		GapOpenPenalty:           10,
		GapExtensionPenalty:      2,
		IndelCodonOpeningBonus:   0,
		IndelCodonExtensionBonus: 2,
		PositionalIndelScores: map[int][2]int{
			8:  [2]int{-6, 0},
			9:  [2]int{-3, 0},
			10: [2]int{6, 0},
			11: [2]int{-3, 0},
			12: [2]int{-6, 0},
		},
		SupportsPositionalIndelScores: true,
	})
	aln, _ = NewAlignment(nseq, ASEQ, handler)
	result = aln.GetReport()
	//          1  2  3  4  5  6  7  8  9 10       11 12 13 14 15 16 17 18 19
	//        ACAGTRTTAGTAGGACCTACACCTttttttGCCAACATAATTGGAAGAAATCTGTTGACYCAG
	expect = "::::::::::::::::::::::::......++++++:::::::::::::::::::::::::::"
	if result.ControlLine != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result.ControlLine)
	}
	nseq = n.ReadString("ACAGTRTTAGTAGGACCTACACCTGCCAACgggATAATTGGAAGAAATCTGTTGACYCAG")
	aln, _ = NewAlignment(nseq, ASEQ, handler)
	result = aln.GetReport()
	//          1  2  3  4  5  6  7  8  9 10    11 12 13 14 15 16 17 18 19
	//        ACAGTRTTAGTAGGACCTACACCTGCCAACgggATAATTGGAAGAAATCTGTTGACYCAG
	expect = "::::::::::::::::::::::::...:::+++:::::::::::::::::::::::::::"
	if result.ControlLine != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result.ControlLine)
	}
	nseq = n.ReadString("ACAGTRTTAGTAGGACCTACACCTGCCAACATArrrATTGGAAGAAATCTGTTGACYCAG")
	aln, _ = NewAlignment(nseq, ASEQ, handler)
	result = aln.GetReport()
	//          1  2  3  4  5  6  7  8  9 10    11 12 13 14 15 16 17 18 19
	//        ACAGTRTTAGTAGGACCTACACCTGCCAACATArrrATTGGAAGAAATCTGTTGACYCAG
	expect = "::::::::::::::::::::::::...:::+++...::::::::::::::::::::::::"
	if result.ControlLine != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result.ControlLine)
	}
	nseq = n.ReadString("ACAGTRTTAGTAGGACCTACACCTGCCAACATAATTGrrrGAAGAAATCTGTTGACYCAG")
	aln, _ = NewAlignment(nseq, ASEQ, handler)
	result = aln.GetReport()
	//          1  2  3  4  5  6  7  8  9 10 11 12 13    14 15 16 17 18 19
	//        ACAGTRTTAGTAGGACCTACACCTGCCAACATAATTGrrrGAAGAAATCTGTTGACYCAG
	expect = "::::::::::::::::::::::::...:::::::::...+++::::::::::::::::::"
	if result.ControlLine != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result.ControlLine)
	}
	nseq = n.ReadString("ACAGTRTTAGTAGGACCcccTACACCTGCCAACATAATTGGAAGAAATCTGTTGACYCAG")
	aln, _ = NewAlignment(nseq, ASEQ, handler)
	result = aln.GetReport()
	//          1  2  3  4  5  6     7  8  9 10 11 12 13 14 15 16 17 18 19
	//        ACAGTRTTAGTAGGACCcccTACACCTGCCAACATAATTGGAAGAAATCTGTTGACYCAG
	expect = "::::::::::::::::::+++::::::...::::::::::::::::::::::::::::::"
	if result.ControlLine != expect {
		t.Errorf(MSG_NOT_EQUAL, expect, result.ControlLine)
	}

}
