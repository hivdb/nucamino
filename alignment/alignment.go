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
	maxScorePos         tPos
	stopCodonPenalty    int
	nSeqLen             int
	aSeqLen             int
	gapOpenPenalty      int
	gapExtensionPenalty int
	directionMatrix     []int
	codonInsPositions   []int
}

func NewAlignment(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, gapOpenPenalty int, gapExtensionPenalty int) *Alignment {
	nSeqLen := len(nSeq)
	aSeqLen := len(aSeq)
	typedPosLen := scoreTypeCount * (nSeqLen + 1) * (aSeqLen + 1)
	codonInsLen := (nSeqLen + 1) * (aSeqLen + 1)
	directionMatrix := make([]int, typedPosLen)
	for i := 0; i < typedPosLen; i++ {
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
	//nSeqLen := self.nSeqLen
	//aSeqLen := self.aSeqLen
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
	endMtIdx := self.getMatrixIndex(GENERAL, self.maxScorePos)
	/*maxScore := negInf
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
	}*/
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

func (self *Alignment) setCachedScore(scoreType tScoreType, pos tPos, prevMatrixIdx int) {
	mtIdx := self.getMatrixIndex(scoreType, pos)
	self.directionMatrix[mtIdx] = prevMatrixIdx
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

func (self *Alignment) calcExtInsScore(pos tPos,
	gScore10 int, iScore10 int) (int, int) {
	score := negInf
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
		pos10 := posShift(pos, -1, 0)
		cand1 := iScore10 + gapPenalty1
		cand2 := gScore10 + gapPenalty2
		if cand1 >= cand2 {
			score, prevMatrixIdx /*, control*/ = cand1, self.getMatrixIndex(EXT_INS, pos10) //, "+"
		} else {
			score, prevMatrixIdx /*, control*/ = cand2, self.getMatrixIndex(GENERAL, pos10) //, "+"
		}
	}
	return score, prevMatrixIdx
}

func (self *Alignment) calcDelScore(pos tPos,
	gScore01 int, gScore11 int, gScore21 int,
	dScore01 int, dScore11 int) (int, int) {
	score := negInf
	q := self.gapOpenPenalty
	r := self.gapExtensionPenalty
	var prevMatrixIdx int
	//var control string
	if pos.n > 0 && pos.a == 0 {
		score, prevMatrixIdx = CalcGapScore(0, q, r), self.getMatrixIndex(GENERAL, tPos{0, 0})
		//control = strings.Repeat("+", pos.n)
	} else if pos.n > 0 {
		curNA := self.getNA(pos.n)
		var prevNA n.NucleicAcid
		curAA := self.getAA(pos.a)
		pos01 := posShift(pos, 0, -1)
		pos11 := posShift(pos, -1, -1)
		pos21 := posShift(pos, -2, -1)
		score = negInf
		//if pos.n < self.nSeqLen {
		cand := dScore01 + CalcGapScore(3, 0, r)
		if cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos01) //, "---"
		}

		cand = gScore01 + CalcGapScore(3, q, r)
		if cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos01) //, "---"
		}

		cand = gScore11 +
			self.CalcMutationScore(curNA, n.N, n.N, curAA) +
			CalcGapScore(2, q, r)
		if cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos11) //, ".--"
		}

		cand = gScore11 +
			self.CalcMutationScore(n.N, curNA, n.N, curAA) +
			CalcGapScore(1, q, r)*2
		if cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos11) //, "-.-"
		}

		if pos.n > 1 {
			prevNA = self.getNA(pos.n - 1)
			cand = dScore11 +
				self.CalcMutationScore(n.N, curNA, n.N, curAA) +
				CalcGapScore(2, q, r)
			if cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos11) //, "-.-"
			}

			cand = gScore21 +
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
	return score, prevMatrixIdx
}

func (self *Alignment) calcScore(
	pos tPos,
	gScore11 int, gScore21 int, gScore31 int,
	iScore00 int,
	dScore00 int, dScore11 int, dScore21 int) (int, int) {
	var prevMatrixIdx int
	score := negInf
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
		if cand := /* #1 */ gScore11 +
			self.CalcMutationScore(n.N, n.N, curNA, curAA) +
			CalcGapScore(2, q, r); cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos11) //, "--."
		}
		if pos.n > 1 {
			prevNA = self.getNA(pos.n - 1)
			if cand := /* #2 */ gScore21 +
				self.CalcMutationScore(prevNA, n.N, curNA, curAA) +
				CalcGapScore(1, q, r); cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos21) //, ".-."
			}
			if cand := /* #3 */ gScore21 +
				self.CalcMutationScore(n.N, prevNA, curNA, curAA) +
				CalcGapScore(1, q, r); cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos21) //, "-.."
			}
			if cand := /* #7 */ dScore11 +
				self.CalcMutationScore(n.N, n.N, curNA, curAA) +
				CalcGapScore(2, 0, r); cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, pos11) //, "--."
			}
		}
		if pos.n > 2 {
			prevNA2 = self.getNA(pos.n - 2)
			if cand := /* #4 */ gScore31 +
				self.CalcMutationScore(prevNA2, prevNA, curNA, curAA); cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, pos31) //, "..."
			}

			if cand := /* #8 */ dScore21 +
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
		if cand := /* #9 */ iScore00; cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(EXT_INS, pos) //, ""
		}
		if cand := /* #10 */ dScore00; cand >= score {
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
	return score, prevMatrixIdx
}

func (self *Alignment) CalcScore() int {
	maxScore := negInf
	maxScorePos := tPos{0, 0}
	gScores := make([]int, self.nSeqLen+1)
	dScores := make([]int, self.nSeqLen+1)
	gScoresCur := make([]int, self.nSeqLen+1)
	dScoresCur := make([]int, self.nSeqLen+1)

	for j := 0; j <= self.aSeqLen; j++ {
		gScore10 := negInf
		iScore10 := negInf
		for i := 0; i <= self.nSeqLen; i++ {
			pos := tPos{i, j}

			gScore01 := gScores[i]
			dScore01 := dScores[i]
			gScore11, gScore21, gScore31 := negInf, negInf, negInf
			dScore11, dScore21 := negInf, negInf
			if i > 0 {
				gScore11 = gScores[i-1]
				dScore11 = dScores[i-1]
			}
			if i > 1 {
				gScore21 = gScores[i-2]
				dScore21 = dScores[i-2]
			}
			if i > 2 {
				gScore31 = gScores[i-3]
			}
			iScore00, iPrevMtIdx := self.calcExtInsScore(pos,
				gScore10, iScore10)

			self.setCachedScore(EXT_INS, pos, iPrevMtIdx)

			dScore00, dPrevMtIdx := self.calcDelScore(pos,
				gScore01, gScore11, gScore21,
				dScore01, dScore11)

			dScoresCur[i] = dScore00
			self.setCachedScore(DEL, pos, dPrevMtIdx)

			gScore00, gPrevMtIdx := self.calcScore(pos,
				gScore11, gScore21, gScore31,
				iScore00,
				dScore00, dScore11, dScore21)

			gScoresCur[i] = gScore00
			self.setCachedScore(GENERAL, pos, gPrevMtIdx)

			if gScore00 >= maxScore {
				maxScore = gScore00
				maxScorePos = pos
			}

			gScore10 = gScore00
			iScore10 = iScore00
		}

		gScores, gScoresCur = gScoresCur, gScores
		dScores, dScoresCur = dScoresCur, dScores

	}
	self.maxScorePos = maxScorePos
	return maxScore
}
