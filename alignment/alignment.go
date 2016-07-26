package alignment

import (
	a "../types/amino"
	n "../types/nucleic"
	"fmt"
	//"github.com/chrislusf/glow/flow"
	"strings"
)

type tScoreType int

const (
	GENERAL   tScoreType = iota
	EXT_INS              // E
	DEL                  // F
	CODON_INS            // C
)

const scoreTypeCount = 4

func scoreTypeToString(scoreType tScoreType) string {
	return map[tScoreType]string{
		GENERAL:   "GENERAL",
		EXT_INS:   "EXT_INS",
		DEL:       "DEL",
		CODON_INS: "CODON_INS",
	}[scoreType]
}

type tPos struct {
	n int
	a int
}

type tTypedPos struct {
	tPos
	scoreType tScoreType
}

type tScoredPos struct {
	tPos
	score int
}

func posShift(pos tPos, nDelta int, aDelta int) tPos {
	return tPos{pos.n + nDelta, pos.a + aDelta}
}

type Alignment struct {
	nSeq                []n.NucleicAcid
	aSeq                []a.AminoAcid
	nSeqLen             int
	aSeqLen             int
	gapOpenPenalty      int
	gapExtensionPenalty int
	scoreMatrix         []int
	directionMatrix     []int
	controlMatrix       []string
	typedPosExist       []bool
	codonInsPositions   []int
	codonInsPosExist    []bool
}

func NewAlignment(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, gapOpenPenalty int, gapExtensionPenalty int) *Alignment {
	nSeqLen := len(nSeq)
	aSeqLen := len(aSeq)
	typedPosLen := scoreTypeCount * (nSeqLen + 1) * (aSeqLen + 1)
	codonInsLen := (nSeqLen + 1) * (aSeqLen + 1)
	return &Alignment{
		nSeq:                nSeq,
		aSeq:                aSeq,
		nSeqLen:             nSeqLen,
		aSeqLen:             aSeqLen,
		gapOpenPenalty:      gapOpenPenalty,
		gapExtensionPenalty: gapExtensionPenalty,
		scoreMatrix:         make([]int, typedPosLen),
		directionMatrix:     make([]int, typedPosLen),
		controlMatrix:       make([]string, typedPosLen),
		typedPosExist:       make([]bool, typedPosLen),
		codonInsPositions:   make([]int, codonInsLen),
		codonInsPosExist:    make([]bool, codonInsLen),
	}
}

const negInf = -int((^uint(0))>>1) - 1

func (self *Alignment) GetCalcTimes() int {
	counter := 0
	for _, cond := range self.typedPosExist {
		if cond {
			counter++
		}
	}
	return counter
}

func (self *Alignment) GetScorePath() []tScoredPos {
	nSeqLen := self.nSeqLen
	aSeqLen := self.aSeqLen
	result := make([]tScoredPos, 100, 100)
	print("          ")
	for i := 0; i <= aSeqLen; i++ {
		fmt.Printf("%4d ", i)
	}
	print("\n               ")
	for _, aa := range self.aSeq {
		fmt.Printf("%4s ", a.ToString(aa))
	}
	print("\n")
	//endPos := self.getMatrixIndex(GENERAL, tPos{nSeqLen, aSeqLen})
	maxScore := negInf
	for i := 0; i <= nSeqLen; i++ {
		fmt.Printf("%4d ", i)
		if i > 0 {
			fmt.Printf("%4s ", n.ToString(self.getNA(i)))
		} else {
			print("     ")
		}
		for j := 0; j <= aSeqLen; j++ {
			tmpPos := self.getMatrixIndex(GENERAL, tPos{i, j})
			score := self.scoreMatrix[tmpPos]
			fmt.Printf("%4d ", int(score))
			if score >= maxScore {
				maxScore = score
				//endPos = tmpPos
			}
		}
		print("\n")
	}
	return result
}

func (self *Alignment) getMatrixIndex(scoreType tScoreType, pos tPos) int {
	return (self.aSeqLen+1)*(pos.n+int(scoreType)*(self.nSeqLen+1)) + pos.a
}

func (self *Alignment) getCachedScore(scoreType tScoreType, pos tPos) (int, bool) {
	mtIdx := self.getMatrixIndex(scoreType, pos)
	score, present := 0, self.typedPosExist[mtIdx]
	if present {
		score = self.scoreMatrix[mtIdx]
	}
	return score, present
}

func (self *Alignment) setCachedScore(scoreType tScoreType, pos tPos, score int, prevMatrixIdx int, ctl string) {
	mtIdx := self.getMatrixIndex(scoreType, pos)
	self.scoreMatrix[mtIdx] = score
	self.directionMatrix[mtIdx] = prevMatrixIdx
	self.controlMatrix[mtIdx] = ctl
	self.typedPosExist[mtIdx] = true
}

func (self *Alignment) getNA(nPos int) n.NucleicAcid {
	return self.nSeq[nPos-1]
}

func (self *Alignment) getAA(aPos int) a.AminoAcid {
	return self.aSeq[aPos-1]
}

func (self *Alignment) calcCodonInsPosition(pos tPos) int {
	posIdx := self.getMatrixIndex(GENERAL, pos) // hack: GENERAL (0) means no offset
	resultPos, present := self.codonInsPositions[posIdx], self.codonInsPosExist[posIdx]
	if present {
		return resultPos
	}
	if pos.n <= 3 && pos.a > 0 {
		resultPos = 0
	}

	r := self.gapExtensionPenalty
	prevPos := posShift(pos, -1, 0)
	codonInsScoreCand1 := self.calcCodonInsScore(pos)
	codonInsScoreCand2 := self.calcCodonInsScore(prevPos) + CalcGapScore(1, 0, r)

	if pos.n >= 4 && pos.a > 0 {
		if codonInsScoreCand1 > codonInsScoreCand2 {
			resultPos = pos.n - 3
		} else if codonInsScoreCand1 == codonInsScoreCand2 {
			resultPos = self.calcCodonInsPosition(prevPos)
		}
	}
	self.codonInsPositions[posIdx] = resultPos
	self.codonInsPosExist[posIdx] = true
	return resultPos
}

func (self *Alignment) calcCodonInsScore(pos tPos) int {
	score, present := self.getCachedScore(CODON_INS, pos)
	if present {
		return score
	}
	q := self.gapOpenPenalty
	r := self.gapExtensionPenalty
	var prevMatrixIdx int
	var control string
	if pos.n <= 3 && pos.a > 0 {
		score, prevMatrixIdx = negInf, self.getMatrixIndex(GENERAL, tPos{0, 0})
		control = strings.Repeat("-", pos.a*3-pos.n) + strings.Repeat(".", pos.n)
	} else if pos.a > 0 {
		prevPos1 := posShift(pos, -1, 0)
		prevPos2 := posShift(pos, -4, -1)
		cand1 := self.calcCodonInsScore(prevPos1) + CalcGapScore(1, 0, r)
		cand2 := self.calcScore(prevPos2) + CalcGapScore(1, q, r)
		if cand1 > cand2 {
			score, prevMatrixIdx, control = cand1, self.getMatrixIndex(CODON_INS, prevPos1), ""
		} else {
			score, prevMatrixIdx, control = cand2, self.getMatrixIndex(GENERAL, prevPos2), ""
		}
	}
	self.setCachedScore(CODON_INS, pos, score, prevMatrixIdx, control)
	return score
}

func (self *Alignment) calcExtInsScore(pos tPos) int {
	score, present := self.getCachedScore(EXT_INS, pos)
	if present {
		return score
	}
	q := self.gapOpenPenalty
	r := self.gapExtensionPenalty
	var prevMatrixIdx int
	var control string
	if pos.n == 0 && pos.a > 0 {
		score, prevMatrixIdx = CalcGapScore(0, q, r), self.getMatrixIndex(GENERAL, tPos{0, 0})
		control = strings.Repeat("---", pos.a)
	} else {
		gapPenalty1 := 0
		gapPenalty2 := 0
		if pos.a > 0 && pos.a < self.aSeqLen {
			gapPenalty1 = CalcGapScore(1, 0, r)
			gapPenalty2 = CalcGapScore(1, q, r)
		}
		prevPos := posShift(pos, -1, 0)
		cand1 := self.calcExtInsScore(prevPos) + gapPenalty1
		cand2 := self.calcScore(prevPos) + gapPenalty2
		if cand1 >= cand2 {
			score, prevMatrixIdx, control = cand1, self.getMatrixIndex(EXT_INS, prevPos), "+"
		} else {
			score, prevMatrixIdx, control = cand2, self.getMatrixIndex(GENERAL, prevPos), "+"
		}
	}
	self.setCachedScore(EXT_INS, pos, score, prevMatrixIdx, control)
	return score
}

func (self *Alignment) calcDelScore(pos tPos) int {
	score, present := self.getCachedScore(DEL, pos)
	if present {
		return score
	}
	q := self.gapOpenPenalty
	r := self.gapExtensionPenalty
	var prevMatrixIdx int
	var control string
	if pos.n > 0 && pos.a == 0 {
		score, prevMatrixIdx = CalcGapScore(0, q, r), self.getMatrixIndex(GENERAL, tPos{0, 0})
		control = strings.Repeat("+", pos.n)
	} else {
		curNA := self.getNA(pos.n)
		var prevNA n.NucleicAcid
		curAA := self.getAA(pos.a)
		pos01 := posShift(pos, 0, -1)
		pos11 := posShift(pos, -1, -1)
		pos21 := posShift(pos, -2, -1)
		score = negInf
		if pos.n < self.nSeqLen {
			cand := self.calcDelScore(pos01) + CalcGapScore(3, 0, r)
			if cand >= score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(DEL, pos01), "---"
			}

			cand = self.calcScore(pos01) + CalcGapScore(3, q, r)
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos01), "---"
			}

			cand = self.calcScore(pos11) +
				CalcMutationScore(curNA, n.N, n.N, curAA) +
				CalcGapScore(2, q, r)
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos11), ".--"
			}

			cand = self.calcScore(pos11) +
				CalcMutationScore(n.N, curNA, n.N, curAA) +
				CalcGapScore(1, q, r)*2
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos11), "-.-"
			}

			if pos.n > 1 {
				prevNA = self.getNA(pos.n - 1)
				cand = self.calcDelScore(pos11) +
					CalcMutationScore(n.N, curNA, n.N, curAA) +
					CalcGapScore(2, q, r)
				if cand >= score {
					score, prevMatrixIdx, control = cand, self.getMatrixIndex(DEL, pos11), "-.-"
				}

				cand = self.calcScore(pos21) +
					CalcMutationScore(prevNA, curNA, n.N, curAA) +
					CalcGapScore(1, q, r)
				if cand >= score {
					score, prevMatrixIdx, control = cand, self.getMatrixIndex(DEL, pos21), "..-"
				}
			}
		} else {
			cand1 := self.calcDelScore(pos01)
			cand2 := self.calcScore(pos01)
			if cand1 >= cand2 {
				score, prevMatrixIdx, control = cand1, self.getMatrixIndex(DEL, pos01), "-"
			} else {
				score, prevMatrixIdx, control = cand2, self.getMatrixIndex(GENERAL, pos01), "-"
			}
		}
	}
	self.setCachedScore(DEL, pos, score, prevMatrixIdx, control)
	return score
}

func (self *Alignment) calcScore(pos tPos) int {
	//fmt.Printf("%s\n", pos)
	score, present := self.getCachedScore(GENERAL, pos)
	if present {
		return score
	}
	var prevMatrixIdx int
	var control string
	if pos.n == 0 || pos.a == 0 {
		score, prevMatrixIdx = 0, self.getMatrixIndex(GENERAL, pos)
	} else {
		q := self.gapOpenPenalty
		r := self.gapExtensionPenalty
		curNA := self.getNA(pos.n)
		var prevNA, prevNA2, prevNA3 n.NucleicAcid
		curAA := self.getAA(pos.a)
		codonInsPos := self.calcCodonInsPosition(pos)
		score = negInf
		pos11 := posShift(pos, -1, -1)
		pos21 := posShift(pos, -2, -1)
		pos31 := posShift(pos, -3, -1)
		pos41 := posShift(pos, -4, -1)
		cand := /* #1 */ self.calcScore(pos11) +
			CalcMutationScore(n.N, n.N, curNA, curAA) +
			CalcGapScore(2, q, r)

		if cand > score {
			score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos11), "--."
		}
		if pos.n > 1 {
			prevNA = self.getNA(pos.n - 1)
			cand = /* #2 */ self.calcScore(pos21) +
				CalcMutationScore(prevNA, n.N, curNA, curAA) +
				CalcGapScore(1, q, r)
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos21), ".-."
			}

			cand = /* #3 */ self.calcScore(pos21) +
				CalcMutationScore(n.N, prevNA, curNA, curAA) +
				CalcGapScore(1, q, r)
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos21), "-.."
			}

			cand = /* #7 */ self.calcDelScore(pos11) +
				CalcMutationScore(n.N, n.N, curNA, curAA) +
				CalcGapScore(2, 0, r)
			if cand >= score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(DEL, pos11), "--."
			}
		}
		if pos.n > 2 {
			prevNA2 = self.getNA(pos.n - 2)
			cand = /* #4 */ self.calcScore(pos31) +
				CalcMutationScore(prevNA2, prevNA, curNA, curAA)
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos31), "..."
			}

			cand = /* #8 */ self.calcDelScore(pos21) +
				CalcMutationScore(n.N, prevNA, curNA, curAA) +
				CalcGapScore(1, 0, r)
			if cand >= score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(DEL, pos21), "-.."
			}
		}
		if pos.n > 3 {
			prevNA3 = self.getNA(pos.n - 3)
			cand = /* #5 */ self.calcScore(pos41) +
				CalcMutationScore(prevNA3, prevNA, curNA, curAA) +
				CalcGapScore(1, q, r)
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos41), ".+.."
			}

			cand = /* #6 */ self.calcScore(pos41) +
				CalcMutationScore(prevNA3, prevNA2, curNA, curAA) +
				CalcGapScore(1, q, r)
			if cand > score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos41), "..+."
			}
		}
		cand = /* #9 */ self.calcExtInsScore(pos)
		if cand >= score {
			score, prevMatrixIdx, control = cand, self.getMatrixIndex(EXT_INS, pos), ""
		}
		cand = /* #10 */ self.calcDelScore(pos)
		if cand >= score {
			score, prevMatrixIdx, control = cand, self.getMatrixIndex(DEL, pos), ""
		}
		if codonInsPos > 0 {
			cand = /* #11 */ self.calcCodonInsScore(pos) +
				CalcMutationScore(
					self.getNA(codonInsPos),
					prevNA, curNA, curAA)
			if cand >= score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(CODON_INS, pos),
					"."+strings.Repeat("+", pos.n-1-codonInsPos)+".."
			}
			cand = /* #12 */ self.calcCodonInsScore(pos) +
				CalcMutationScore(
					self.getNA(codonInsPos),
					self.getNA(codonInsPos+1),
					curNA, curAA)
			if cand >= score {
				score, prevMatrixIdx, control = cand, self.getMatrixIndex(CODON_INS, pos),
					".."+strings.Repeat("+", pos.n-1-codonInsPos)+"."
			}
		}
	}
	self.setCachedScore(GENERAL, pos, score, prevMatrixIdx, control)
	return score
}

func (self *Alignment) CalcScore() int {
	return self.calcScore(tPos{self.nSeqLen, self.aSeqLen})
}
