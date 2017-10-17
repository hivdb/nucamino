package alignment

import (
	"errors"
	// "fmt"
	h "github.com/hivdb/nucamino/scorehandler/general"
	a "github.com/hivdb/nucamino/types/amino"
	f "github.com/hivdb/nucamino/types/frameshift"
	m "github.com/hivdb/nucamino/types/mutation"
	n "github.com/hivdb/nucamino/types/nucleic"
	sortutil "github.com/pmylund/sortutil"
	// "os"
	"strings"
)

type tScoreType int

const (
	GENERAL tScoreType = iota
	INS                // E
	DEL                // F
	//CODON_INS            // C
)

const negInf = -int((^uint(0))>>1) - 1
const scoreTypeCount = 3

func (self tScoreType) ToString() string {
	return map[tScoreType]string{
		GENERAL: "GENERAL",
		INS:     "INS",
		DEL:     "DEL",
		//CODON_INS: "CODON_INS",
	}[self]
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
	scoreHandler                  h.GeneralScoreHandler
	scoreMatrix                   []int
	q                             int
	r                             int
	supportPositionalIndel        bool
	constIndelCodonOpeningScore   int
	constIndelCodonExtensionScore int
	boundaryOnly                  bool
	isSimpleAlignment             bool
}

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

func NewAlignment(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, scoreHandler h.GeneralScoreHandler) (*Alignment, error) {
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
		scoreMatrix:                   make([]int, 0),
		supportPositionalIndel:        supportPositionalIndel,
		constIndelCodonOpeningScore:   constIndelCodonOpeningScore,
		constIndelCodonExtensionScore: constIndelCodonExtensionScore,
	}

	success := result.align()
	if success {
		return result, nil
	} else {
		return nil, errors.New("sequence misaligned")
	}
}

func (self *Alignment) GetCalcInfo() (int, int) {
	counter := 0
	for idx, mtIdx := range self.scoreMatrix {
		if idx%2 == 0 {
			// score
			continue
		}
		if mtIdx != -1 {
			counter++
		}
	}
	return counter, len(self.scoreMatrix) / 2
}

func padRightSpace(str string, length int) string {
	return str + strings.Repeat(" ", (length-len(str)))
}

func (self *Alignment) GetReport() AlignmentReport {
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
		endMtIdx                         = self.getMatrixIndex(GENERAL, self.endPosN, self.endPosA)
		startMtIdx                       = 0
		singleMtLen                      = (self.nSeqLen + 1) * (self.aSeqLen + 1) * 2
	)

	if self.isSimpleAlignment {
		startMtIdx = self.getMatrixIndex(GENERAL, 0, 0)
	} else {
		var (
			score     int
			curMtIdx  = endMtIdx
			prevMtIdx = -1
			prevScore = self.maxScore
		)
		for curMtIdx >= 0 && curMtIdx != prevMtIdx {
			score = self.scoreMatrix[curMtIdx]
			if score <= prevScore {
				prevScore = score
				startMtIdx = curMtIdx
			}
			prevMtIdx = curMtIdx
			curMtIdx = self.scoreMatrix[curMtIdx+1]
		}
	}
	// fmt.Fprintf(os.Stderr, "1,1 - %d,%d (%d)\n", self.endPosA, self.endPosN, singleMtLen)

	for endMtIdx%singleMtLen >= startMtIdx%singleMtLen {
		// score := self.scoreMatrix[endMtIdx]
		scoreType, posN, posA := self.getTypedPos(endMtIdx)
		// fmt.Fprintf(os.Stderr, "%s (%d,%d,%d,%d)\n", lastScoreType, LastPosN, posN, LastPosA, posA)
		if lastAA == 0 && lastNA == 0 {
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
				partialALine += padRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
				partialCLine += ":::"
			} else {
				if mutation != nil {
					partialCLine += mutation.Control
					if mutation.IsDeletion {
						partialNLine += "   "
						partialALine += padRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
					} else {
						partialNLine += mutation.CodonText
						partialALine += padRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
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
		if self.isSimpleAlignment {
			endMtIdx = self.getMatrixIndex(GENERAL, posN-3, posA-1)
		} else {
			endMtIdx = self.scoreMatrix[endMtIdx+1]
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
	/*print(aLine, ")\n")
	print(cLine, ")\n")
	print(nLine, ")\n")
	print("First AA: ", firstAA, "\n")
	print("First NA: ", firstNA, "\n")
	print("Last AA: ", lastAA, "\n")
	print("Last NA: ", lastNA, "\n")
	for _, mutation := range mutList {
		print(mutation.ToString(), ", ")
	}
	print("\n")
	for _, frameshift := range fsList {
		print(frameshift.ToString(), ", ")
	}
	print("\n")
	print("Total Score: ", self.maxScore, "\n")*/
	sortutil.Reverse(mutList)
	sortutil.Reverse(fsList)
	sortutil.Reverse(siteList)
	return AlignmentReport{
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
}

func (self *Alignment) getMatrixIndex(scoreType tScoreType, posN int, posA int) int {
	return 2 * ((self.aSeqLen+1)*(posN+int(scoreType)*(self.nSeqLen+1)) + posA)
}

func (self *Alignment) getTypedPos(matrixIndex int) (scoreType tScoreType, posN int, posA int) {
	matrixIndex /= 2
	posA = matrixIndex % (self.aSeqLen + 1)
	nTotal := matrixIndex / (self.aSeqLen + 1)
	posN = nTotal % (self.nSeqLen + 1)
	scoreType = tScoreType(nTotal / (self.nSeqLen + 1))
	return
}

func (self *Alignment) setCachedScore(scoreType tScoreType, posN int, posA int, score int, prevMatrixIdx int) {
	mtIdx := self.getMatrixIndex(scoreType, posN, posA)
	self.scoreMatrix[mtIdx] = score
	self.scoreMatrix[mtIdx+1] = prevMatrixIdx
}

func (self *Alignment) getNA(nPos int) n.NucleicAcid {
	return self.nSeq[nPos-1]
}

func (self *Alignment) getAA(aPos int) a.AminoAcid {
	return self.aSeq[aPos-1]
}

func (self *Alignment) align() bool {
	var endPosN, endPosA, simplesCount int
	self.boundaryOnly = true
	// set boundary for nSeq, so we don't have to build a huge nSeqLen * aSeqLen matrix
	endPosN, endPosA, self.maxScore, simplesCount = self.calcScoreMainForward()
	self.nSeq = self.nSeq[:endPosN]
	self.nSeqLen = len(self.nSeq)
	self.aSeq = self.aSeq[:endPosA]
	self.aSeqLen = len(self.aSeq)
	startPosN, startPosA, _ := self.calcScoreMainBackward()
	if endPosN <= startPosN {
		return false
	}
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
		return true
	}
	typedPosLen := scoreTypeCount * (self.nSeqLen + 1) * (self.aSeqLen + 1) * 2
	self.scoreMatrix = make([]int, typedPosLen)
	self.endPosN, self.endPosA, self.maxScore, _ = self.calcScoreMainForward()
	return true
}
