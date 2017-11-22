package alignment

import (
	n "github.com/hivdb/nucamino/types/nucleic"
)

func (self *Alignment) calcExtInsScoreBackward(
	posN int, posA int,
	gScore30 int, iScore30 int,
	gScore20 int, gScore10 int) int {
	var (
		insOpeningScore   int
		insExtensionScore int
		cand              int
		score             = negInf
		sh                = self.scoreHandler
		q                 = self.q
		r                 = self.r
	)
	//var control string
	if posN == self.nSeqLen && posA < self.aSeqLen {
		score = 0
		//control = strings.Repeat("---", pos.a)
	} else {
		score = negInf
		if posA == 1 {
			q, r, insOpeningScore, insExtensionScore = 0, 0, 0, 0
		} else {
			if self.supportPositionalIndel {
				insOpeningScore, insExtensionScore = sh.GetPositionalIndelCodonScore(posA-1, true)
			} else {
				insOpeningScore = self.constIndelCodonOpeningScore
				insExtensionScore = self.constIndelCodonExtensionScore
			}
		}
		if posN < self.nSeqLen-3 {
			if cand = iScore30 + r + r + r + insExtensionScore; cand > score {
				score = cand // "+++"
			}
			if cand = gScore30 + q + r + r + r + insOpeningScore + insExtensionScore; cand > score {
				score = cand // "+++"
			}
		}
		if posN < self.nSeqLen-2 {
			if cand = gScore20 + q + r + r; cand > score {
				score = cand // "++"
			}
		}
		if cand = gScore10 + q + r; cand > score {
			score = cand // "+"
		}
	}
	return score
}

func (self *Alignment) calcDelScoreBackward(
	posN int, posA int,
	gScore01 int, gScore11 int, gScore21 int,
	dScore01 int, dScore11 int) int {
	var (
		score = negInf
		sh    = self.scoreHandler
	)
	//var control string
	if posN < self.nSeqLen && posA == self.aSeqLen {
		score = 0
		//control = strings.Repeat("+", pos.n)
	} else if posN < self.nSeqLen {
		var (
			delOpeningScore   int
			delExtensionScore int
			q, q2, r, r2      = self.q, self.q, self.r, self.r
		)
		score = negInf
		if posN == 1 {
			q, r, delOpeningScore, delExtensionScore = 0, 0, 0, 0
		} else {
			if self.supportPositionalIndel {
				delOpeningScore, delExtensionScore = sh.GetPositionalIndelCodonScore(posA-1, false)
			} else {
				delOpeningScore = self.constIndelCodonOpeningScore
				delExtensionScore = self.constIndelCodonExtensionScore
			}
		}
		if cand := dScore01 + r + r + r + delExtensionScore; cand >= score {
			score = cand // "---"
		}

		if cand := gScore01 + q + r + r + r + delOpeningScore + delExtensionScore; cand > score {
			score = cand // "---"
		}

		if cand := gScore11 + q + r + r; cand > score {
			score = cand // "--."
		}

		if cand := gScore11 + q + r + q2 + r2; cand > score {
			score = cand // "-.-"
		}

		if posN < self.nSeqLen-1 {
			if cand := dScore11 + q + r + r2; cand >= score {
				score = cand // "-.-"
			}

			if cand := gScore21 + q + r; cand >= score {
				score = cand // "-.."
			}
		}
	}
	return score
}

func (self *Alignment) calcScoreBackward(
	posN int, posA int,
	gScore11 int, gScore21 int, gScore31 int,
	iScore00 int,
	dScore00 int, dScore11 int, dScore21 int) int {
	score := negInf
	if posN == self.nSeqLen || posA == self.aSeqLen {
		score = 0
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
			score = cand // ".--"
		}
		if posN < self.nSeqLen-1 {
			prevNA = self.getNA(posN + 1)
			if cand := /* #2 */ gScore21 + q + r; cand > score {
				score = cand // ".-."
			}
			if cand := /* #3 */ gScore21 + q + r; cand > score {
				score = cand // "..-"
			}
			if cand := /* #7 */ dScore11 + r + r; cand >= score {
				score = cand // ".--"
			}
		}
		if posN < self.nSeqLen-2 {
			prevNA2 = self.getNA(posN + 2)
			tmpScore, present := sh.GetCachedSubstitutionScore(posA, curNA, prevNA, prevNA2, curAA)
			if !present {
				tmpScore = sh.GetSubstitutionScoreNoCache(posA, curNA, prevNA, prevNA2, curAA)
			}
			if cand := /* #4 */ gScore31 + tmpScore; cand > score {
				score = cand // "..."
			}

			if cand := /* #8 */ dScore21 + r; cand >= score {
				score = cand // "..-"
			}
		}
		if cand := /* #9 */ iScore00; cand >= score {
			score = cand // ""
		}
		if cand := /* #10 */ dScore00; cand >= score {
			score = cand // ""
		}
	}
	return score
}

func (self *Alignment) calcScoreMainBackward() (int, int, int) {
	var (
		maxScore     = negInf
		maxScorePosN = self.nSeqLen
		maxScorePosA = self.aSeqLen
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
	)

	for j := self.aSeqLen; j >= 1; j-- {
		gScore30, iScore30 = negInf, negInf
		gScore20, iScore20 = negInf, negInf
		gScore10, iScore10 = negInf, negInf
		for i := self.nSeqLen; i >= 1; i-- {
			gScore01 = gScores[i]
			dScore01 = dScores[i]
			gScore11, gScore21, gScore31 = negInf, negInf, negInf
			dScore11, dScore21 = negInf, negInf
			if i < self.nSeqLen {
				gScore11 = gScores[i+1]
				dScore11 = dScores[i+1]
			}
			if i < self.nSeqLen-1 {
				gScore21 = gScores[i+2]
				dScore21 = dScores[i+2]
			}
			if i < self.nSeqLen-2 {
				gScore31 = gScores[i+3]
			}
			iScore00 = self.calcExtInsScoreBackward(
				i, j, gScore30,
				iScore30, gScore20, gScore10)

			dScore00 = self.calcDelScoreBackward(
				i, j,
				gScore01, gScore11, gScore21,
				dScore01, dScore11)

			dScoresCur[i] = dScore00

			gScore00 = self.calcScoreBackward(
				i, j,
				gScore11, gScore21, gScore31,
				iScore00,
				dScore00, dScore11, dScore21)

			gScoresCur[i] = gScore00

			if (i == 1 || j == 1) && gScore00 > maxScore {
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
