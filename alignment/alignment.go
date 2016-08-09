package alignment

import (
	s "../scorehandler"
	a "../types/amino"
	f "../types/frameshift"
	m "../types/mutation"
	n "../types/nucleic"
	"strings"
	//"fmt"
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
	nSeq                   []n.NucleicAcid
	aSeq                   []a.AminoAcid
	nSeqLen                int
	aSeqLen                int
	maxScorePosN           int
	maxScorePosA           int
	maxScore               int
	scoreHandler           s.ScoreHandler
	scoreMatrix            []int
	q                      int
	r                      int
	supportPositionalIndel bool
	//codonInsPositions   []int
}

type AlignmentReport struct {
	FirstAA     int
	FirstNA     int
	LastAA      int
	LastNA      int
	Mutations   []m.Mutation
	FrameShifts []f.FrameShift
}

func NewAlignment(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, scoreHandler s.ScoreHandler) *Alignment {
	nSeqLen := len(nSeq)
	aSeqLen := len(aSeq)
	typedPosLen := scoreTypeCount * (nSeqLen + 1) * (aSeqLen + 1) * 2
	//codonInsLen := (nSeqLen + 1) * (aSeqLen + 1)
	scoreMatrix := make([]int, typedPosLen)
	//codonInsPositions := make([]int, codonInsLen)
	//for i := range codonInsPositions {
	//	codonInsPositions[i] = -1
	//}
	result := &Alignment{
		nSeq:         nSeq,
		aSeq:         aSeq,
		nSeqLen:      nSeqLen,
		aSeqLen:      aSeqLen,
		scoreHandler: scoreHandler,
		scoreMatrix:  scoreMatrix,
		q:            scoreHandler.GetGapOpeningScore(),
		r:            scoreHandler.GetGapExtensionScore(),
		supportPositionalIndel: scoreHandler.IsPositionalIndelScoreSupported(),
		//codonInsPositions:   codonInsPositions,
	}
	result.calcScoreMain()
	return result
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
	//nSeqLen := self.nSeqLen
	//aSeqLen := self.aSeqLen
	//print("          ")
	//for i := 0; i <= aSeqLen; i++ {
	//	fmt.Printf("%4d ", i)
	//}
	//print("\n               ")
	//for _, aa := range self.aSeq {
	//	fmt.Printf("%4s ", a.ToString(aa))
	//}
	//print("\n")
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
	var (
		nLine, aLine, cLine              string
		firstAA, lastAA, firstNA, lastNA int
		mutList                          = make([]m.Mutation, 0, 10)
		fsList                           = make([]f.FrameShift, 0, 3)
		LastPosN                         = -1
		LastPosA                         = -1
		lastScoreType                    = GENERAL
		hasUnprocessedNAs                = false
		endMtIdx                         = self.getMatrixIndex(GENERAL, self.maxScorePosN, self.maxScorePosA)
		prevMtIdx                        = -1
		curMtIdx                         = endMtIdx
		prevScore                        = endMtIdx
		startMtIdx                       = 0
		singleMtLen                      = (self.nSeqLen + 1) * (self.aSeqLen + 1) * 2
	)

	for curMtIdx >= 0 && curMtIdx != prevMtIdx {
		score := self.scoreMatrix[curMtIdx]
		if score <= prevScore {
			prevScore = score
			startMtIdx = curMtIdx
		}
		prevMtIdx = curMtIdx
		curMtIdx = self.scoreMatrix[curMtIdx+1]
	}

	for endMtIdx%singleMtLen >= startMtIdx%singleMtLen {
		//score := self.scoreMatrix[endMtIdx]
		scoreType, posN, posA := self.getTypedPos(endMtIdx)
		//fmt.Printf("%s (%d,%d,%d,%d,%d)\n", lastScoreType, LastPosN, posN, LastPosA, posA, score)
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
				mutation = m.MakeMutation(posA+1, self.nSeq[posN:LastPosN], self.aSeq[posA])
				frameshift = f.MakeFrameShift(posA+1, self.nSeq[posN:LastPosN])
				if mutation != nil {
					mutList = append(mutList, *mutation)
				}
				if frameshift != nil {
					fsList = append(fsList, *frameshift)
				}
			}
			/* those are only for generate three lines */
			if LastPosA > posA && LastPosN-posN > 2 && mutation == nil {
				partialNLine += n.WriteString(self.nSeq[posN : posN+3])
				partialALine += padRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
				partialCLine += ":::"
			} else {
				if mutation != nil {
					partialCLine += mutation.GetControl()
					if mutation.IsDeletion() {
						partialNLine += "   "
						partialALine += padRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
					} else {
						codon := mutation.GetCodon()
						partialNLine += codon.ToString()
						partialALine += padRightSpace(a.WriteString(self.aSeq[posA:LastPosA]), 3)
						if mutation.IsInsertion() {
							for _, insCodon := range mutation.GetInsertedCodons() {
								partialNLine += insCodon.ToString()
								partialALine += "   "
							}
						}
					}
				}
			}
			if frameshift != nil && frameshift.IsInsertion() {
				partialNLine += n.WriteString(frameshift.GetNucleicAcids())
				partialALine += strings.Repeat(" ", frameshift.GetGapLength())
				partialCLine += strings.Repeat("+", frameshift.GetGapLength())
			}
			/* end */

			nLine = partialNLine + nLine
			aLine = partialALine + aLine
			cLine = partialCLine + cLine
		}
		//control = self.controlMatrix[endMtIdx]
		endMtIdx = self.scoreMatrix[endMtIdx+1]
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
	print(aLine, ")\n")
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
	print("Total Score: ", self.maxScore, "\n")
	return AlignmentReport{
		FirstAA:     firstAA,
		FirstNA:     firstNA,
		LastAA:      lastAA,
		LastNA:      lastNA,
		Mutations:   mutList,
		FrameShifts: fsList,
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

func (self *Alignment) calcExtInsScore(
	posN int, posA int,
	gScore30 int, iScore30 int,
	gScore10 int, iScore10 int) (int, int) {
	var (
		cand, prevMatrixIdx int
		score               = negInf
		sh                  = self.scoreHandler
		q                   = self.q
		r                   = self.r
	)
	//var control string
	if posN == 0 && posA > 0 {
		score, prevMatrixIdx = q, self.getMatrixIndex(GENERAL, 0, 0)
		//control = strings.Repeat("---", pos.a)
	} else {
		//gapPenalty1 := 0
		//gapPenalty2 := 0
		//if pos.a > 0 && pos.a < self.aSeqLen {
		//}
		score = negInf
		if self.supportPositionalIndel && posN > 3 {
			insScore := sh.GetIndelScore(posA, true)
			if cand = iScore30 + r + r + r + insScore; cand > score {
				score, prevMatrixIdx = cand, self.getMatrixIndex(INS, posN-3, posA)
			}
			if cand = gScore30 + q + r + r + r + insScore; cand > score {
				score, prevMatrixIdx = cand, self.getMatrixIndex(GENERAL, posN-3, posA)
			}
		}
		if cand = iScore10 + r; cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(INS, posN-1, posA) //, "+"
		}
		if cand = gScore10 + q + r; cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-1, posA) //, "+"
		}
	}
	return score, prevMatrixIdx
}

func (self *Alignment) calcDelScore(
	posN int, posA int,
	gScore01 int, gScore11 int, gScore21 int,
	dScore01 int, dScore11 int) (int, int) {
	var (
		prevMatrixIdx int
		score         = negInf
		sh            = self.scoreHandler
	)
	//var control string
	if posN > 0 && posA == 0 {
		score, prevMatrixIdx = self.q, self.getMatrixIndex(GENERAL, 0, 0)
		//control = strings.Repeat("+", pos.n)
	} else if posN > 0 {
		var (
			prevNA      n.NucleicAcid
			curNA       = self.getNA(posN)
			curAA       = self.getAA(posA)
			q           = self.q
			r           = self.r
			mutScoreN0N = sh.GetSubstitutionScore(posA, n.N, curNA, n.N, curAA)
			delScore    = 0
		)
		score = negInf
		//if pos.n < self.nSeqLen {
		if self.supportPositionalIndel {
			delScore = sh.GetIndelScore(posA, false)
		}
		if cand := dScore01 + r + r + r + delScore; cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, posN, posA-1) //, "---"
		}

		if cand := gScore01 + q + r + r + r + delScore; cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN, posA-1) //, "---"
		}

		if cand := gScore11 +
			sh.GetSubstitutionScore(posA, curNA, n.N, n.N, curAA) + q + r + r; cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-1, posA-1) //, ".--"
		}

		if cand := gScore11 + mutScoreN0N + q + r + q + r; cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-1, posA-1) //, "-.-"
		}

		if posN > 1 {
			prevNA = self.getNA(posN - 1)
			if cand := dScore11 + mutScoreN0N + q + r + r; cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, posN-1, posA-1) //, "-.-"
			}

			if cand := gScore21 +
				sh.GetSubstitutionScore(posA, prevNA, curNA, n.N, curAA) + q + r; cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-2, posA-1) //, "..-"
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
	posN int, posA int,
	gScore11 int, gScore21 int, gScore31 int,
	iScore00 int,
	dScore00 int, dScore11 int, dScore21 int) (int, int) {
	var (
		prevMatrixIdx int
		score         = negInf
	)
	if posN == 0 || posA == 0 {
		score, prevMatrixIdx = 0, self.getMatrixIndex(GENERAL, posN, posA)
	} else {
		var (
			prevNA, prevNA2/*, prevNA3*/ n.NucleicAcid
			mutScoreN10 int
			sh          = self.scoreHandler
			q           = self.q
			r           = self.r
			curNA       = self.getNA(posN)
			curAA       = self.getAA(posA)
			mutScoreNN0 = sh.GetSubstitutionScore(posA, n.N, n.N, curNA, curAA)
		)
		//codonInsPos := self.calcCodonInsPosition(pos)
		score = negInf
		//pos41 := posShift(pos, -4, -1)
		if cand := /* #1 */ gScore11 + mutScoreNN0 + q + r + r; cand > score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-1, posA-1) //, "--."
		}
		if posN > 1 {
			prevNA = self.getNA(posN - 1)
			mutScoreN10 = sh.GetSubstitutionScore(posA, n.N, prevNA, curNA, curAA)
			if cand := /* #2 */ gScore21 +
				sh.GetSubstitutionScore(posA, prevNA, n.N, curNA, curAA) + q + r; cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-2, posA-1) //, ".-."
			}
			if cand := /* #3 */ gScore21 + mutScoreN10 + q + r; cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-2, posA-1) //, "-.."
			}
			if cand := /* #7 */ dScore11 + mutScoreNN0 + r + r; cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, posN-1, posA-1) //, "--."
			}
		}
		if posN > 2 {
			prevNA2 = self.getNA(posN - 2)
			if cand := /* #4 */ gScore31 +
				sh.GetSubstitutionScore(posA, prevNA2, prevNA, curNA, curAA); cand > score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(GENERAL, posN-3, posA-1) //, "..."
			}

			if cand := /* #8 */ dScore21 + mutScoreN10 + r; cand >= score {
				score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, posN-2, posA-1) //, "-.."
			}
		}
		if cand := /* #9 */ iScore00; cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(INS, posN, posA) //, ""
		}
		if cand := /* #10 */ dScore00; cand >= score {
			score, prevMatrixIdx /*, control*/ = cand, self.getMatrixIndex(DEL, posN, posA) //, ""
		}
	}
	return score, prevMatrixIdx
}

func (self *Alignment) calcScoreMain() {
	var (
		maxScore     = negInf
		maxScorePosN = 0
		maxScorePosA = 0
		gScores      = make([]int, self.nSeqLen+1)
		dScores      = make([]int, self.nSeqLen+1)
		gScoresCur   = make([]int, self.nSeqLen+1)
		dScoresCur   = make([]int, self.nSeqLen+1)

		gScore30, gScore20 int
		gScore00, gScore10 int
		gScore01, gScore11 int
		gScore21, gScore31 int

		iScore30, iScore20 int
		iScore00, iScore10 int

		dScore00, dScore01 int
		dScore11, dScore21 int
		prevMtIdx          int
	)

	for j := 0; j <= self.aSeqLen; j++ {
		gScore30, iScore30 = negInf, negInf
		gScore20, iScore20 = negInf, negInf
		gScore10, iScore10 = negInf, negInf
		for i := 0; i <= self.nSeqLen; i++ {

			gScore01 = gScores[i]
			dScore01 = dScores[i]
			gScore11, gScore21, gScore31 = negInf, negInf, negInf
			dScore11, dScore21 = negInf, negInf
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
			iScore00, prevMtIdx = self.calcExtInsScore(
				i, j, gScore30, iScore30, gScore10, iScore10)

			self.setCachedScore(INS, i, j, iScore00, prevMtIdx)

			dScore00, prevMtIdx = self.calcDelScore(
				i, j,
				gScore01, gScore11, gScore21,
				dScore01, dScore11)

			dScoresCur[i] = dScore00
			self.setCachedScore(DEL, i, j, dScore00, prevMtIdx)

			gScore00, prevMtIdx = self.calcScore(
				i, j,
				gScore11, gScore21, gScore31,
				iScore00,
				dScore00, dScore11, dScore21)

			gScoresCur[i] = gScore00
			self.setCachedScore(GENERAL, i, j, gScore00, prevMtIdx)

			if gScore00 >= maxScore {
				maxScore = gScore00
				maxScorePosN = i
				maxScorePosA = j
			}

			gScore30, gScore20, gScore10 = gScore20, gScore10, gScore00
			iScore30, iScore20, iScore10 = iScore20, iScore10, iScore00
		}

		gScores, gScoresCur = gScoresCur, gScores
		dScores, dScoresCur = dScoresCur, dScores

	}
	self.maxScorePosN = maxScorePosN
	self.maxScorePosA = maxScorePosA
	self.maxScore = maxScore
}
