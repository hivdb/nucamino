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
	)
	//var control string
	if posN == 0 && posA > 0 {
		score = q
		if !self.boundaryOnly {
			prevMatrixIdx = self.getMatrixIndex(GENERAL, 0, 0)
		}
		//control = strings.Repeat("---", pos.a)
	} else {
		score = negInf
		if self.supportPositionalIndel {
			insOpeningScore, insExtensionScore = sh.GetPositionalIndelCodonScore(posA, true)
		} else {
			insOpeningScore = self.constIndelCodonOpeningScore
			insExtensionScore = self.constIndelCodonExtensionScore
		}
		if posN > 3 {
			if cand = iScore30 + r + r + r + insExtensionScore; cand > score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(INS, posN-3, posA) //, "+++"
				}
			}
			if cand = gScore30 + q + r + r + r + insOpeningScore + insExtensionScore; cand > score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-3, posA) //, "+++"
				}
			}
		}
		if posN > 2 {
			if cand = gScore20 + q + r + r; cand > score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-2, posA) //, "++"
				}
			}
		}
		if cand = gScore10 + q + r; cand > score {
			score = cand
			if !self.boundaryOnly {
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
	)
	//var control string
	if posN > 0 && posA == 0 {
		score = self.q
		if !self.boundaryOnly {
			prevMatrixIdx = self.getMatrixIndex(GENERAL, 0, 0)
		}
		//control = strings.Repeat("+", pos.n)
	} else if posN > 0 {
		var (
			delOpeningScore   int
			delExtensionScore int
			prevNA            n.NucleicAcid
			curNA             = self.getNA(posN)
			curAA             = self.getAA(posA)
			q                 = self.q
			r                 = self.r
			mutScoreN0N       = sh.GetSubstitutionScore(posA, n.N, curNA, n.N, curAA)
		)
		score = negInf
		if self.supportPositionalIndel {
			delOpeningScore, delExtensionScore = sh.GetPositionalIndelCodonScore(posA, false)
		} else {
			delOpeningScore = self.constIndelCodonOpeningScore
			delExtensionScore = self.constIndelCodonExtensionScore
		}
		if cand := dScore01 + r + r + r + delExtensionScore; cand >= score {
			score = cand
			if !self.boundaryOnly {
				prevMatrixIdx = self.getMatrixIndex(DEL, posN, posA-1) //, "---"
			}
		}

		if cand := gScore01 + q + r + r + r + delOpeningScore + delExtensionScore; cand > score {
			score = cand
			if !self.boundaryOnly {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN, posA-1) //, "---"
			}
		}

		if cand := gScore11 +
			sh.GetSubstitutionScore(posA, curNA, n.N, n.N, curAA) + q + r + r; cand > score {
			score = cand
			if !self.boundaryOnly {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-1, posA-1) //, ".--"
			}
		}

		if cand := gScore11 + mutScoreN0N + q + r + q + r; cand > score {
			score = cand
			if !self.boundaryOnly {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-1, posA-1) //, "-.-"
			}
		}

		if posN > 1 {
			prevNA = self.getNA(posN - 1)
			if cand := dScore11 + mutScoreN0N + q + r + r; cand >= score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(DEL, posN-1, posA-1) //, "-.-"
				}
			}

			if cand := gScore21 +
				sh.GetSubstitutionScore(posA, prevNA, curNA, n.N, curAA) + q + r; cand >= score {
				score = cand
				if !self.boundaryOnly {
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
	dScore00 int, dScore11 int, dScore21 int) (int, int) {
	var (
		prevMatrixIdx int
		score         = negInf
	)
	if posN == 0 || posA == 0 {
		score = 0
		if !self.boundaryOnly {
			prevMatrixIdx = self.getMatrixIndex(GENERAL, posN, posA)
		}
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
		score = negInf
		if cand := /* #1 */ gScore11 + mutScoreNN0 + q + r + r; cand > score {
			score = cand
			if !self.boundaryOnly {
				prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-1, posA-1) //, "--."
			}
		}
		if posN > 1 {
			prevNA = self.getNA(posN - 1)
			mutScoreN10 = sh.GetSubstitutionScore(posA, n.N, prevNA, curNA, curAA)
			if cand := /* #2 */ gScore21 +
				sh.GetSubstitutionScore(posA, prevNA, n.N, curNA, curAA) + q + r; cand > score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-2, posA-1) //, ".-."
				}
			}
			if cand := /* #3 */ gScore21 + mutScoreN10 + q + r; cand > score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-2, posA-1) //, "-.."
				}
			}
			if cand := /* #7 */ dScore11 + mutScoreNN0 + r + r; cand >= score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(DEL, posN-1, posA-1) //, "--."
				}
			}
		}
		if posN > 2 {
			prevNA2 = self.getNA(posN - 2)
			if cand := /* #4 */ gScore31 +
				sh.GetSubstitutionScore(posA, prevNA2, prevNA, curNA, curAA); cand > score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(GENERAL, posN-3, posA-1) //, "..."
				}
			}

			if cand := /* #8 */ dScore21 + mutScoreN10 + r; cand >= score {
				score = cand
				if !self.boundaryOnly {
					prevMatrixIdx = self.getMatrixIndex(DEL, posN-2, posA-1) //, "-.."
				}
			}
		}
		if cand := /* #9 */ iScore00; cand >= score {
			score = cand
			if !self.boundaryOnly {
				prevMatrixIdx = self.getMatrixIndex(INS, posN, posA) //, ""
			}
		}
		if cand := /* #10 */ dScore00; cand >= score {
			score = cand
			if !self.boundaryOnly {
				prevMatrixIdx = self.getMatrixIndex(DEL, posN, posA) //, ""
			}
		}
	}
	return score, prevMatrixIdx
}

func (self *Alignment) calcScoreMainForward() (int, int, int) {
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
			iScore00, prevMtIdx = self.calcExtInsScoreForward(
				i, j, gScore30,
				iScore30, gScore20, gScore10)

			if !self.boundaryOnly {
				self.setCachedScore(INS, i, j, iScore00, prevMtIdx)
			}

			dScore00, prevMtIdx = self.calcDelScoreForward(
				i, j,
				gScore01, gScore11, gScore21,
				dScore01, dScore11)

			dScoresCur[i] = dScore00

			if !self.boundaryOnly {
				self.setCachedScore(DEL, i, j, dScore00, prevMtIdx)
			}

			gScore00, prevMtIdx = self.calcScoreForward(
				i, j,
				gScore11, gScore21, gScore31,
				iScore00,
				dScore00, dScore11, dScore21)

			gScoresCur[i] = gScore00

			if !self.boundaryOnly {
				self.setCachedScore(GENERAL, i, j, gScore00, prevMtIdx)
			}

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
	return maxScorePosN, maxScorePosA, maxScore
}
