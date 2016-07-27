package alignment

import (
	a "../types/amino"
	c "../types/codon"
	n "../types/nucleic"
	"fmt"
	//"strings"
)

type tScoreType int

const (
	GENERAL tScoreType = iota
	EXT_INS            // E
	DEL                // F
	//CODON_INS            // C
)

const negInf = -int((^uint(0))>>1) - 1
const scoreTypeCount = 3

func scoreTypeToString(scoreType tScoreType) string {
	return map[tScoreType]string{
		GENERAL: "GENERAL",
		EXT_INS: "EXT_INS",
		DEL:     "DEL",
		//CODON_INS: "CODON_INS",
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
	stopCodonPenalty    int
	nSeqLen             int
	aSeqLen             int
	gapOpenPenalty      int
	gapExtensionPenalty int
	scoreMatrix         []int
	directionMatrix     []int
	codonInsPositions   []int
}

func NewAlignment(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, gapOpenPenalty int, gapExtensionPenalty int) *Alignment {
	nSeqLen := len(nSeq)
	aSeqLen := len(aSeq)
	typedPosLen := scoreTypeCount * (nSeqLen + 1) * (aSeqLen + 1)
	codonInsLen := (nSeqLen + 1) * (aSeqLen + 1)
	scoreMatrix := make([]int, typedPosLen)
	directionMatrix := make([]int, typedPosLen)
	for i := 0; i < typedPosLen; i++ {
		scoreMatrix[i] = negInf
		directionMatrix[i] = -1
	}
	codonInsPositions := make([]int, codonInsLen)
	for i := range codonInsPositions {
		codonInsPositions[i] = -1
	}
	return &Alignment{
		nSeq:                nSeq,
		aSeq:                aSeq,
		stopCodonPenalty:    5,
		nSeqLen:             nSeqLen,
		aSeqLen:             aSeqLen,
		gapOpenPenalty:      gapOpenPenalty,
		gapExtensionPenalty: gapExtensionPenalty,
		scoreMatrix:         scoreMatrix,
		directionMatrix:     directionMatrix,
		codonInsPositions:   codonInsPositions,
	}
}

func (self *Alignment) GetCalcInfo() (int, int) {
	counter := 0
	for _, prevIdx := range self.directionMatrix {
		if prevIdx != -1 {
			counter++
		}
	}
	return counter, len(self.directionMatrix)
}

func (self *Alignment) GetScorePath() []tScoredPos {
	nSeqLen := self.nSeqLen
	aSeqLen := self.aSeqLen
	result := make([]tScoredPos, 100, 100)
	//print("          ")
	//for i := 0; i <= aSeqLen; i++ {
	//	fmt.Printf("%4d ", i)
	//}
	//print("\n               ")
	//for _, aa := range self.aSeq {
	//	fmt.Printf("%4s ", a.ToString(aa))
	//}
	//print("\n")
	endMtIdx := self.getMatrixIndex(GENERAL, tPos{0, 0})
	maxScore := negInf
	for i := 0; i <= nSeqLen; i++ {
		//fmt.Printf("%4d ", i)
		//if i > 0 {
		//	fmt.Printf("%4s ", n.ToString(self.getNA(i)))
		//} else {
		//	print("     ")
		//}
		for j := 0; j <= aSeqLen; j++ {
			tmpMtIdx := self.getMatrixIndex(GENERAL, tPos{i, j})
			score := self.scoreMatrix[tmpMtIdx]
			//fmt.Printf("%4d ", int(score))
			if score >= maxScore {
				maxScore = score
				endMtIdx = tmpMtIdx
			}
		}
		//print("\n")
	}
	lastN, lastA := -1, -1
	var nLine, aLine, cLine string
	for endMtIdx >= 0 {
		//score := self.scoreMatrix[endMtIdx]
		_, pos := self.getTypedPos(endMtIdx)
		//fmt.Printf("%s (%d,%d) %d\n", scoreType, pos.n, pos.a, score)
		if lastN > -1 && lastA > -1 {
			padding := (lastA - pos.a) * 3
			if padding < lastN-pos.n {
				padding = lastN - pos.n
			}
			format := fmt.Sprintf("%%%ds", padding)
			nLine = fmt.Sprintf(format, n.WriteString(self.nSeq[pos.n:lastN])) + nLine
			aLine = fmt.Sprintf(format, a.WriteString(self.aSeq[pos.a:lastA])) + aLine
			cLine = c.GetFinalControlLine(self.nSeq[pos.n:lastN], self.aSeq[pos.a:lastA]) + cLine
		}
		//control = self.controlMatrix[endMtIdx]
		endMtIdx = self.directionMatrix[endMtIdx]
		lastN = pos.n
		lastA = pos.a
		if pos.n == 0 || pos.a == 0 {
			break
		}
	}
	print(aLine, ")\n")
	print(cLine, ")\n")
	print(nLine, ")\n")
	return result
}

func (self *Alignment) getMatrixIndex(scoreType tScoreType, pos tPos) int {
	return (self.aSeqLen+1)*(pos.n+int(scoreType)*(self.nSeqLen+1)) + pos.a
}

func (self *Alignment) getTypedPos(matrixIndex int) (tScoreType, tPos) {
	a := matrixIndex % (self.aSeqLen + 1)
	nTotal := matrixIndex / (self.aSeqLen + 1)
	n := nTotal % (self.nSeqLen + 1)
	scoreType := tScoreType(nTotal / (self.nSeqLen + 1))
	return scoreType, tPos{n, a}
}

func (self *Alignment) getCachedScore(scoreType tScoreType, pos tPos) (int, bool) {
	mtIdx := self.getMatrixIndex(scoreType, pos)
	score := self.scoreMatrix[mtIdx]
	present := score != negInf
	return score, present
}

func (self *Alignment) setCachedScore(scoreType tScoreType, pos tPos, score int, prevMatrixIdx int /*, ctl string*/) {
	mtIdx := self.getMatrixIndex(scoreType, pos)
	self.scoreMatrix[mtIdx] = score
	self.directionMatrix[mtIdx] = prevMatrixIdx
	//self.controlMatrix[mtIdx] = ctl
}

func (self *Alignment) getNA(nPos int) n.NucleicAcid {
	return self.nSeq[nPos-1]
}

func (self *Alignment) getAA(aPos int) a.AminoAcid {
	return self.aSeq[aPos-1]
}

/*func (self *Alignment) calcCodonInsPosition(pos tPos) int {
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
}*/

/*func (self *Alignment) calcCodonInsScore(pos tPos) int {
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
}*/

func (self *Alignment) calcExtInsScore(pos tPos) int {
	score, present := self.getCachedScore(EXT_INS, pos)
	if present {
		return score
	}
	q := self.gapOpenPenalty
	r := self.gapExtensionPenalty
	var prevMatrixIdx int
	//var control string
	if pos.n == 0 && pos.a > 0 {
		score, prevMatrixIdx = CalcGapScore(0, q, r), self.getMatrixIndex(GENERAL, tPos{0, 0})
		//control = strings.Repeat("---", pos.a)
	} else {
		//gapPenalty1 := 0
		//gapPenalty2 := 0
		//if pos.a > 0 && pos.a < self.aSeqLen {
		gapPenalty1 := CalcGapScore(1, 0, r)
		gapPenalty2 := CalcGapScore(1, q, r)
		//}
		prevPos := posShift(pos, -1, 0)
		cand1 := self.calcExtInsScore(prevPos) + gapPenalty1
		cand2 := self.calcScore(prevPos) + gapPenalty2
		if cand1 >= cand2 {
			score, prevMatrixIdx /*, control*/ = cand1, self.getMatrixIndex(EXT_INS, prevPos) //, "+"
		} else {
			score, prevMatrixIdx /*, control*/ = cand2, self.getMatrixIndex(GENERAL, prevPos) //, "+"
		}
	}
	self.setCachedScore(EXT_INS, pos, score, prevMatrixIdx /*, control*/)
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
	//var control string
	if pos.n > 0 && pos.a == 0 {
		score, prevMatrixIdx = CalcGapScore(0, q, r), self.getMatrixIndex(GENERAL, tPos{0, 0})
		//control = strings.Repeat("+", pos.n)
	} else {
		curNA := self.getNA(pos.n)
		var prevNA n.NucleicAcid
		curAA := self.getAA(pos.a)
		pos01 := posShift(pos, 0, -1)
		pos11 := posShift(pos, -1, -1)
		pos21 := posShift(pos, -2, -1)
		score = negInf
		//if pos.n < self.nSeqLen {
		cand := self.calcDelScore(pos01) + CalcGapScore(3, 0, r)
		if cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos01) //, "---"
		}

		cand = self.calcScore(pos01) + CalcGapScore(3, q, r)
		if cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos01) //, "---"
		}

		cand = self.calcScore(pos11) +
			self.CalcMutationScore(curNA, n.N, n.N, curAA) +
			CalcGapScore(2, q, r)
		if cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos11) //, ".--"
		}

		cand = self.calcScore(pos11) +
			self.CalcMutationScore(n.N, curNA, n.N, curAA) +
			CalcGapScore(1, q, r)*2
		if cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos11) //, "-.-"
		}

		if pos.n > 1 {
			prevNA = self.getNA(pos.n - 1)
			cand = self.calcDelScore(pos11) +
				self.CalcMutationScore(n.N, curNA, n.N, curAA) +
				CalcGapScore(2, q, r)
			if cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos11) //, "-.-"
			}

			cand = self.calcScore(pos21) +
				self.CalcMutationScore(prevNA, curNA, n.N, curAA) +
				CalcGapScore(1, q, r)
			if cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos21) //, "..-"
			}
		}
		/*} else {
			// pos.n == self.nSeqLen
			cand1 := self.calcDelScore(pos01)
			cand2 := self.calcScore(pos01)
			if cand1 >= cand2 {
				score, prevMatrixIdx, control = cand1, self.getMatrixIndex(DEL, pos01), "---"
			} else {
				score, prevMatrixIdx, control = cand2, self.getMatrixIndex(GENERAL, pos01), "---"
			}
		}*/
	}
	self.setCachedScore(DEL, pos, score, prevMatrixIdx /*, control*/)
	return score
}

func (self *Alignment) calcScore(pos tPos) int {
	//fmt.Printf("%s\n", pos)
	score, present := self.getCachedScore(GENERAL, pos)
	if present {
		return score
	}
	var prevMatrixIdx int
	//var control string
	if pos.n == 0 || pos.a == 0 {
		score, prevMatrixIdx = 0, self.getMatrixIndex(GENERAL, pos)
	} else {
		q := self.gapOpenPenalty
		r := self.gapExtensionPenalty
		curNA := self.getNA(pos.n)
		var prevNA, prevNA2 /*, prevNA3*/ n.NucleicAcid
		curAA := self.getAA(pos.a)
		//codonInsPos := self.calcCodonInsPosition(pos)
		score = negInf
		pos11 := posShift(pos, -1, -1)
		pos21 := posShift(pos, -2, -1)
		pos31 := posShift(pos, -3, -1)
		//pos41 := posShift(pos, -4, -1)
		if cand := /* #1 */ self.calcScore(pos11) +
			self.CalcMutationScore(n.N, n.N, curNA, curAA) +
			CalcGapScore(2, q, r); cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos11) //, "--."
		}
		if pos.n > 1 {
			prevNA = self.getNA(pos.n - 1)
			if cand := /* #2 */ self.calcScore(pos21) +
				self.CalcMutationScore(prevNA, n.N, curNA, curAA) +
				CalcGapScore(1, q, r); cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos21) //, ".-."
			}
			if cand := /* #3 */ self.calcScore(pos21) +
				self.CalcMutationScore(n.N, prevNA, curNA, curAA) +
				CalcGapScore(1, q, r); cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos21) //, "-.."
			}
			if cand := /* #7 */ self.calcDelScore(pos11) +
				self.CalcMutationScore(n.N, n.N, curNA, curAA) +
				CalcGapScore(2, 0, r); cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos11) //, "--."
			}
		}
		if pos.n > 2 {
			prevNA2 = self.getNA(pos.n - 2)
			if cand := /* #4 */ self.calcScore(pos31) +
				self.CalcMutationScore(prevNA2, prevNA, curNA, curAA); cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos31) //, "..."
			}

			if cand := /* #8 */ self.calcDelScore(pos21) +
				self.CalcMutationScore(n.N, prevNA, curNA, curAA) +
				CalcGapScore(1, 0, r); cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos21) //, "-.."
			}
		}
		/*
			// This part of code is not necessary for the purpose to align known virus sequence to known reference
			if pos.n > 3 {
				prevNA3 = self.getNA(pos.n - 3)
				#5 cand = self.calcScore(pos41) +
					self.CalcMutationScore(prevNA3, prevNA, curNA, curAA) +
					CalcGapScore(1, q, r)
				if cand > score {
					score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos41), ".+.."
				}

				#6 cand = self.calcScore(pos41) +
					self.CalcMutationScore(prevNA3, prevNA2, curNA, curAA) +
					CalcGapScore(1, q, r)
				if cand > score {
					score, prevMatrixIdx, control = cand, self.getMatrixIndex(GENERAL, pos41), "..+."
				}
			}
		*/
		if cand := /* #9 */ self.calcExtInsScore(pos); cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(EXT_INS, pos) //, ""
		}
		if cand := /* #10 */ self.calcDelScore(pos); cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos) //, ""
		}
		/*
			// This part of code is not necessary for the purpose to align known virus sequence to known reference
			if codonInsPos > 0 {
				#11 cand = self.calcCodonInsScore(pos) +
					self.CalcMutationScore(
						self.getNA(codonInsPos),
						prevNA, curNA, curAA)
				if cand >= score {
					score, prevMatrixIdx, control = cand, self.getMatrixIndex(CODON_INS, pos),
						"."+strings.Repeat("+", pos.n-1-codonInsPos)+".."
				}
				#12 cand = self.calcCodonInsScore(pos) +
					self.CalcMutationScore(
						self.getNA(codonInsPos),
						self.getNA(codonInsPos+1),
						curNA, curAA)
				if cand >= score {
					score, prevMatrixIdx, control = cand, self.getMatrixIndex(CODON_INS, pos),
						".."+strings.Repeat("+", pos.n-1-codonInsPos)+"."
				}
			}
		*/
	}
	self.setCachedScore(GENERAL, pos, score, prevMatrixIdx /*, control*/)
	return score
}

func (self *Alignment) CalcScore() int {
	/*score := negInf
	step := 10

	for i, j := step*3, 0; i <= self.nSeqLen && j <= self.aSeqLen; {
		newScore := self.calcScore(tPos{i, j})

	}*/
	return self.calcScore(tPos{self.nSeqLen, self.aSeqLen})
}
