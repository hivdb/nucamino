package alignment

import (
	n "github.com/hivdb/nucamino/types/nucleic"
)

func (self *Alignment) calcExtInsScoreForward(
	posN int, posA int,
	gScore30 int, iScore30 int,
	gScore20 int, gScore10 int) (int, int) {
	var (
		insOpeningScore     int
		insExtensionScore   int
		cand, prevMatrixIdx int
		score               = negInf
		sh                  = self.scoreHandler
		q                   = self.q
		r                   = self.r
		calcMtIdx           = !self.boundaryOnly
	)
	//var control string
	if posN == 0 && posA > 0 {
		score = 0 // no penalty for initial gaps
		if calcMtIdx {
			prevMatrixIdx = self.getMatrixIndex(GENERAL, 0, 0)
		}
		//control = strings.Repeat("---", pos.a)
	} else {
		score = negInf
		if posA == self.aSeqLen {
			// no penalty for trailing gaps
			r, q, insOpeningScore, insExtensionScore = 0, 0, 0, 0
		} else {
			if self.supportPositionalIndel {
				insOpeningScore, insExtensionScore = sh.GetPositionalIndelCodonScore(posA+self.aSeqOffset, true)
			} else {
				insOpeningScore = self.constIndelCodonOpeningScore
				insExtensionScore = self.constIndelCodonExtensionScore
			}
		}
		if posN > 3 {
			if cand = iScore30 + r + r + r + insExtensionScore; cand > score {
				score = cand
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(INS, posN-3, posA) //, "+++"
				}
			}
			if cand = gScore30 + q + r + r + r + insOpeningScore + insExtensionScore; cand > score {
				score = cand
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-3, posA) //, "+++"
				}
			}
		}
		if posN > 2 {
			if cand = gScore20 + q + r + r; cand > score {
				score = cand
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-2, posA) //, "++"
				}
			}
		}
		if cand = gScore10 + q + r; cand > score {
			score = cand
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-1, posA) //, "+"
			}
		}
	}
	return score, prevMatrixIdx
}

func (self *Alignment) calcDelScoreForward(
	posN int, posA int,
	gScore01 int, gScore11 int, gScore21 int,
	dScore01 int, dScore11 int) (int, int) {
	var (
		prevMatrixIdx int
		score         = negInf
		sh            = self.scoreHandler
		calcMtIdx     = !self.boundaryOnly
	)
	//var control string
	if posN > 0 && posA == 0 {
		score = 0 // no penalty for initial gaps
		if calcMtIdx {
			prevMatrixIdx = self.getMatrixIndex(GENERAL, 0, 0)
		}
		//control = strings.Repeat("+", pos.n)
	} else if posN > 0 {
		var (
			delOpeningScore   int
			delExtensionScore int
			q, q2, r, r2      = self.q, self.q, self.r, self.r
		)
		score = negInf
		if posN == self.nSeqLen {
			// no penalty for trailing gaps
			q, r, delOpeningScore, delExtensionScore = 0, 0, 0, 0
		} else {
			if self.supportPositionalIndel {
				delOpeningScore, delExtensionScore = sh.GetPositionalIndelCodonScore(posA+self.aSeqOffset, false)
			} else {
				delOpeningScore = self.constIndelCodonOpeningScore
				delExtensionScore = self.constIndelCodonExtensionScore
			}
		}
		if cand := dScore01 + r + r + r + delExtensionScore; cand >= score {
			score = cand
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(DEL, posN, posA-1) //, "---"
			}
		}

		if cand := gScore01 + q + r + r + r + delOpeningScore + delExtensionScore; cand > score {
			score = cand
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN, posA-1) //, "---"
			}
		}

		if cand := gScore11 + q + r + r; cand > score {
			score = cand
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-1, posA-1) //, ".--"
			}
		}

		if cand := gScore11 + q2 + r2 + q + r; cand > score {
			score = cand
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-1, posA-1) //, "-.-"
			}
		}

		if posN > 1 {
			if cand := dScore11 + r2 + q + r; cand >= score {
				score = cand
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(DEL, posN-1, posA-1) //, "-.-"
				}
			}

			if cand := gScore21 + q + r; cand >= score {
				score = cand
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-2, posA-1) //, "..-"
				}
			}
		}
	}
	return score, prevMatrixIdx
}

func (self *Alignment) calcScoreForward(
	posN int, posA int,
	gScore11 int, gScore21 int, gScore31 int,
	iScore00 int,
	dScore00 int, dScore11 int, dScore21 int) (int, int, bool) {
	var (
		prevMatrixIdx int
		isSimple      bool
		score         = negInf
		calcMtIdx     = !self.boundaryOnly
	)
	if posN == 0 || posA == 0 {
		score = 0
		if calcMtIdx {
			prevMatrixIdx = self.getMatrixIndex(GENERAL, posN, posA)
		}
	} else {
		var (
			prevNA, prevNA2/*, prevNA3*/ n.NucleicAcid
			sh    = self.scoreHandler
			q     = self.q
			r     = self.r
			curNA = self.getNA(posN)
			curAA = self.getAA(posA)
		)
		score = negInf
		if cand := /* #1 */ gScore11 + q + r + r; cand > score {
			score = cand
			isSimple = false
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-1, posA-1) //, "--."
			}
		}
		if posN > 1 {
			prevNA = self.getNA(posN - 1)
			if cand := /* #2 */ gScore21 + q + r; cand > score {
				score = cand
				isSimple = false
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-2, posA-1) //, ".-."
				}
			}
			if cand := /* #3 */ gScore21 + q + r; cand > score {
				score = cand
				isSimple = false
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-2, posA-1) //, "-.."
				}
			}
			if cand := /* #7 */ dScore11 + r + r; cand >= score {
				score = cand
				isSimple = false
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(DEL, posN-1, posA-1) //, "--."
				}
			}
		}
		if posN > 2 {
			prevNA2 = self.getNA(posN - 2)
			tmpScore, present := sh.GetCachedSubstitutionScore(posA+self.aSeqOffset, prevNA2, prevNA, curNA, curAA)
			if !present {
				tmpScore = sh.GetSubstitutionScoreNoCache(posA+self.aSeqOffset, prevNA2, prevNA, curNA, curAA)
			}
			if cand := /* #4 */ gScore31 + tmpScore; cand > score {
				score = cand
				isSimple = true
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-3, posA-1) //, "..."
				}
			}

			if cand := /* #8 */ dScore21 + r; cand >= score {
				score = cand
				isSimple = false
				if calcMtIdx {
					prevMatrixIdx = self.getMatrixIndex(DEL, posN-2, posA-1) //, "-.."
				}
			}
		}
		if cand := /* #9 */ iScore00; cand >= score {
			score = cand
			isSimple = false
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(INS, posN, posA) //, ""
			}
		}
		if cand := /* #10 */ dScore00; cand >= score {
			score = cand
			isSimple = false
			if calcMtIdx {
				prevMatrixIdx = self.getMatrixIndex(DEL, posN, posA) //, ""
			}
		}
	}
	return score, prevMatrixIdx, isSimple
}

func (self *Alignment) calcScoreMainForward() (int, int, int, int) {
	var (
		maxScore               = negInf
		maxScorePosN           = 0
		maxScorePosA           = 0
		simplesCountAtMaxScore = 0
		gScores                = make([]int, self.nSeqLen+1)
		dScores                = make([]int, self.nSeqLen+1)
		simplesCountMt         = make([]int, self.nSeqLen+1)
		gScoresCur             = make([]int, self.nSeqLen+1)
		dScoresCur             = make([]int, self.nSeqLen+1)
		simplesCountMtCur      = make([]int, self.nSeqLen+1)

		gScore30, gScore20 int
		gScore00, gScore10 int
		gScore01, gScore11 int
		gScore21, gScore31 int
		simplesCount       int

		iScore30, iScore20 int
		iScore00, iScore10 int

		dScore00, dScore01 int
		dScore11, dScore21 int
		prevMtIdx          int
		isSimple           bool
		calcMtIdx          = !self.boundaryOnly
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
				simplesCount = simplesCountMt[i-3]
			}
			iScore00, prevMtIdx = self.calcExtInsScoreForward(
				i, j, gScore30,
				iScore30, gScore20, gScore10)

			if calcMtIdx {
				self.setPrevMatrixIndex(INS, i, j, prevMtIdx)
			}

			dScore00, prevMtIdx = self.calcDelScoreForward(
				i, j,
				gScore01, gScore11, gScore21,
				dScore01, dScore11)

			dScoresCur[i] = dScore00

			if calcMtIdx {
				self.setPrevMatrixIndex(DEL, i, j, prevMtIdx)
			}

			gScore00, prevMtIdx, isSimple = self.calcScoreForward(
				i, j,
				gScore11, gScore21, gScore31,
				iScore00,
				dScore00, dScore11, dScore21)

			if isSimple {
				simplesCount++
				simplesCountMtCur[i] = simplesCount
			} else {
				// reset simplesCount
				simplesCount = 0
				simplesCountMtCur[i] = 0
			}
			gScoresCur[i] = gScore00

			if calcMtIdx {
				self.setPrevMatrixIndex(GENERAL, i, j, prevMtIdx)
			}

			if (i == self.nSeqLen || j == self.aSeqLen) && gScore00 > maxScore {
				maxScore = gScore00
				maxScorePosN = i
				maxScorePosA = j
				simplesCountAtMaxScore = simplesCount
			}

			gScore30, gScore20, gScore10 = gScore20, gScore10, gScore00
			iScore30, iScore20, iScore10 = iScore20, iScore10, iScore00
		}

		gScores, gScoresCur = gScoresCur, gScores
		simplesCountMt, simplesCountMtCur = simplesCountMtCur, simplesCountMt
		dScores, dScoresCur = dScoresCur, dScores

	}
	return maxScorePosN, maxScorePosA, maxScore, simplesCountAtMaxScore
}
