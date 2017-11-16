package alignment

import (
	"errors"
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	f "github.com/hivdb/nucamino/types/frameshift"
	m "github.com/hivdb/nucamino/types/mutation"
	n "github.com/hivdb/nucamino/types/nucleic"
	u "github.com/hivdb/nucamino/utils"
	sortutil "github.com/pmylund/sortutil"
	"strings"
)

type tScoreType int

const (
	GENERAL tScoreType = iota
	INS                // E
	DEL                // F
)

const negInf = -int((^uint(0))>>1) - 1
const scoreTypeCount = 3

type AlignedSite struct {
	PosAA    int
	PosNA    int
	LengthNA int
}

type AlignmentReport struct {
	FirstAA           int
	FirstNA           int
	LastAA            int
	LastNA            int
	Mutations         []m.Mutation
	FrameShifts       []f.FrameShift
	AlignedSites      []AlignedSite
	AminoAcidsLine    string
	ControlLine       string
	NucleicAcidsLine  string
	IsSimpleAlignment bool
}

type Alignment struct {
	nSeq                          []n.NucleicAcid
	aSeq                          []a.AminoAcid
	nSeqLen                       int
	aSeqLen                       int
	nSeqOffset                    int
	aSeqOffset                    int
	endPosN                       int
	endPosA                       int
	maxScore                      int
	scoreHandler                  *h.GeneralScoreHandler
	nwMatrix                      []int
	q                             int
	r                             int
	supportPositionalIndel        bool
	constIndelCodonOpeningScore   int
	constIndelCodonExtensionScore int
	boundaryOnly                  bool
	isSimpleAlignment             bool
	report                        *AlignmentReport
}

func NewAlignment(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, scoreHandler *h.GeneralScoreHandler) (*Alignment, error) {
	nSeqLen := len(nSeq)
	aSeqLen := len(aSeq)
	supportPositionalIndel := scoreHandler.IsPositionalIndelScoreSupported()
	constIndelCodonOpeningScore, constIndelCodonExtensionScore :=
		scoreHandler.GetConstantIndelCodonScore()
	result := &Alignment{
		q:                             scoreHandler.GetGapOpeningScore(),
		r:                             scoreHandler.GetGapExtensionScore(),
		nSeq:                          nSeq,
		aSeq:                          aSeq,
		nSeqLen:                       nSeqLen,
		aSeqLen:                       aSeqLen,
		scoreHandler:                  scoreHandler,
		nwMatrix:                      make([]int, 0),
		supportPositionalIndel:        supportPositionalIndel,
		constIndelCodonOpeningScore:   constIndelCodonOpeningScore,
		constIndelCodonExtensionScore: constIndelCodonExtensionScore,
	}
	ok := result.align()
	if !ok {
		return nil, errors.New("sequence misaligned")
	}
	return result, nil
}

func (self *Alignment) GetReport() *AlignmentReport {
	return self.report
}

func (self *Alignment) generateReport() bool {
	var (
		nLine, aLine, cLine              string
		firstAA, lastAA, firstNA, lastNA int
		mutList                          = make([]m.Mutation, 0, 10)
		fsList                           = make([]f.FrameShift, 0, 3)
		siteList                         = make([]AlignedSite, 0, 50)
		LastPosN                         = -1
		LastPosA                         = -1
		lastScoreType                    = GENERAL
		hasUnprocessedNAs                = false
		startMtIdx                       = self.getMatrixIndex(GENERAL, 0, 0)
		endMtIdx                         = self.getMatrixIndex(GENERAL, self.endPosN, self.endPosA)
		singleMtLen                      = (self.nSeqLen + 1) * (self.aSeqLen + 1) * 2
	)
	for endMtIdx%singleMtLen >= startMtIdx%singleMtLen {
		scoreType, posN, posA := self.getTypedPos(endMtIdx)
		if scoreType != GENERAL && (posN == 0 || posA == 0) {
			break
		}
		if self.isSimpleAlignment {
			endMtIdx = self.getMatrixIndex(GENERAL, posN-3, posA-1)
		} else {
			endMtIdx = self.nwMatrix[endMtIdx]
		}
		if lastAA == 0 && lastNA == 0 {
			if scoreType != GENERAL {
				continue
			}
			nextScoreType, _, _ := self.getTypedPos(endMtIdx)
			if nextScoreType != GENERAL {
				continue
			}
			lastAA, lastNA = posA, posN
		}
		firstAA, firstNA = posA+1, posN+1
		if lastScoreType != INS && LastPosN > -1 && LastPosA > -1 {
			var (
				partialNLine string
				partialALine string
				partialCLine string
				mutation     *m.Mutation
				frameshift   *f.FrameShift
			)

			if LastPosA > posA {
				hasUnprocessedNAs = false
				absPosA := posA + 1 + self.aSeqOffset
				absPosN := posN + 1 + self.nSeqOffset
				lenNA := 3
				mutation = m.MakeMutation(
					absPosA, absPosN,
					self.nSeq[posN:LastPosN], self.aSeq[posA])
				frameshift = f.MakeFrameShift(
					absPosA, absPosN,
					self.nSeq[posN:LastPosN])
				if mutation != nil {
					mutList = append(mutList, *mutation)
					if mutation.IsDeletion {
						lenNA -= 3
					} else if mutation.IsInsertion {
						lenNA += len(mutation.InsertedCodonsText)
					}
				}
				if frameshift != nil {
					fsList = append(fsList, *frameshift)
					if frameshift.IsInsertion {
						lenNA += frameshift.GapLength
					} else {
						lenNA -= frameshift.GapLength
					}
				}
				siteList = append(siteList, AlignedSite{
					PosAA:    absPosA,
					PosNA:    absPosN,
					LengthNA: lenNA,
				})
			}
			/* those are only for generate three lines */
			if LastPosA > posA && LastPosN-posN > 2 && mutation == nil {
				partialNLine += n.WriteString(self.nSeq[posN : posN+3])
				partialALine += u.PadRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
				partialCLine += ":::"
			} else {
				if mutation != nil {
					partialCLine += mutation.Control
					if mutation.IsDeletion {
						partialNLine += "   "
						partialALine += u.PadRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
					} else {
						partialNLine += mutation.CodonText
						partialALine += u.PadRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
						if mutation.IsInsertion {
							for _, insCodon := range mutation.GetInsertedCodons() {
								partialNLine += insCodon.ToString()
								partialALine += "   "
							}
						}
					}
				}
			}
			if frameshift != nil && frameshift.IsInsertion {
				partialNLine += frameshift.NucleicAcidsText
				partialALine += strings.Repeat(" ", frameshift.GapLength)
				partialCLine += strings.Repeat("+", frameshift.GapLength)
			}
			/* end */

			nLine = partialNLine + nLine
			aLine = partialALine + aLine
			cLine = partialCLine + cLine
		}
		if lastScoreType == INS {
			hasUnprocessedNAs = true
		} else if !hasUnprocessedNAs {
			LastPosN = posN
			LastPosA = posA
		}
		lastScoreType = scoreType
		if posN == 0 || posA == 0 {
			break
		}
	}
	sortutil.Reverse(mutList)
	sortutil.Reverse(fsList)
	sortutil.Reverse(siteList)
	self.report = &AlignmentReport{
		FirstAA:           firstAA + self.aSeqOffset,
		FirstNA:           firstNA + self.nSeqOffset,
		LastAA:            lastAA + self.aSeqOffset,
		LastNA:            lastNA + self.nSeqOffset,
		Mutations:         mutList,
		FrameShifts:       fsList,
		AlignedSites:      siteList,
		AminoAcidsLine:    aLine,
		ControlLine:       cLine,
		NucleicAcidsLine:  nLine,
		IsSimpleAlignment: self.isSimpleAlignment,
	}
	return true
}

func (self *Alignment) getMatrixIndex(scoreType tScoreType, posN int, posA int) int {
	return ((self.aSeqLen+1)*(posN+int(scoreType)*(self.nSeqLen+1)) + posA)
}

func (self *Alignment) getTypedPos(matrixIndex int) (scoreType tScoreType, posN int, posA int) {
	posA = matrixIndex % (self.aSeqLen + 1)
	nTotal := matrixIndex / (self.aSeqLen + 1)
	posN = nTotal % (self.nSeqLen + 1)
	scoreType = tScoreType(nTotal / (self.nSeqLen + 1))
	return
}

func (self *Alignment) setPrevMatrixIndex(scoreType tScoreType, posN int, posA int, prevMatrixIdx int) {
	mtIdx := self.getMatrixIndex(scoreType, posN, posA)
	self.nwMatrix[mtIdx] = prevMatrixIdx
}

func (self *Alignment) getNA(nPos int) n.NucleicAcid {
	return self.nSeq[nPos-1]
}

func (self *Alignment) getAA(aPos int) a.AminoAcid {
	return self.aSeq[aPos-1]
}

func (self *Alignment) align() bool {
	var (
		startPosN, startPosA           int
		endPosN, endPosA, simplesCount int
	)
	self.boundaryOnly = true
	// set boundary for nSeq, so we don't have to build a huge nSeqLen * aSeqLen matrix
	endPosN, endPosA, self.maxScore, simplesCount = self.calcScoreMainForward()
	self.nSeq = self.nSeq[:endPosN]
	self.nSeqLen = len(self.nSeq)
	self.aSeq = self.aSeq[:endPosA]
	self.aSeqLen = len(self.aSeq)
	if self.nSeqLen == 0 || self.aSeqLen == 0 {
		return false
	}
	startPosN, startPosA, _ = self.calcScoreMainBackward()
	self.nSeqOffset = startPosN - 1
	self.aSeqOffset = startPosA - 1
	self.boundaryOnly = false
	self.nSeq = self.nSeq[startPosN-1:]
	self.nSeqLen = len(self.nSeq)
	self.aSeq = self.aSeq[startPosA-1:]
	self.aSeqLen = len(self.aSeq)
	if endPosA-self.aSeqOffset == simplesCount {
		self.isSimpleAlignment = true
		self.endPosN = endPosN - self.nSeqOffset
		self.endPosA = endPosA - self.aSeqOffset
		return self.generateReport()
	}
	typedPosLen := scoreTypeCount * (self.nSeqLen + 1) * (self.aSeqLen + 1)
	self.nwMatrix = make([]int, typedPosLen)
	self.endPosN, self.endPosA, self.maxScore, _ = self.calcScoreMainForward()
	return self.generateReport()
}
