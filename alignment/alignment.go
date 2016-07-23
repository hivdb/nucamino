package alignment

import (
	a "../types/amino"
	n "../types/nucleic"
	"fmt"
	//"github.com/chrislusf/glow/flow"
	"math"
)

type tScoreType int

const (
	GENERAL   tScoreType = iota
	EXT_INS              // E
	DEL                  // F
	CODON_INS            // C
)

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
	score float64
}

func posShift(pos tPos, nDelta int, aDelta int) tPos {
	return tPos{pos.n + nDelta, pos.a + aDelta}
}

func flowMax(funcChan chan func() (float64, tTypedPos)) (float64, tTypedPos) {
	score := math.Inf(-1)
	prevPos := tTypedPos{tPos{0, 0}, GENERAL}
	for cb := range funcChan {
		tmpScore, tmpPrevPos := cb()
		if tmpScore > score {
			score = tmpScore
			prevPos = tmpPrevPos
		} else if tmpScore == score &&
			prevPos.scoreType == GENERAL &&
			prevPos.scoreType != tmpPrevPos.scoreType {
			prevPos = tmpPrevPos
		}
	}
	return score, prevPos
}

type Alignment struct {
	nSeq                []n.NucleicAcid
	aSeq                []a.AminoAcid
	gapOpenPenalty      int
	gapExtensionPenalty int
	scoreMatrix         map[tTypedPos]float64
	directionMatrix     map[tTypedPos]tTypedPos
}

func NewAlignment(nSeq []n.NucleicAcid, aSeq []a.AminoAcid, gapOpenPenalty int, gapExtensionPenalty int) *Alignment {
	return &Alignment{
		nSeq:                nSeq,
		aSeq:                aSeq,
		gapOpenPenalty:      gapOpenPenalty,
		gapExtensionPenalty: gapExtensionPenalty,
		scoreMatrix:         make(map[tTypedPos]float64),
		directionMatrix:     make(map[tTypedPos]tTypedPos),
	}
}

func (self *Alignment) GetScorePath() []tScoredPos {
	nSeqLen := self.nSeqLen()
	aSeqLen := self.aSeqLen()
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
	endPos := tTypedPos{tPos{nSeqLen, aSeqLen}, GENERAL}
	maxScore := math.Inf(-1)
	for i := 0; i <= nSeqLen; i++ {
		fmt.Printf("%4d ", i)
		if i > 0 {
			fmt.Printf("%4s ", n.ToString(self.getNA(i)))
		} else {
			print("     ")
		}
		for j := 0; j <= aSeqLen; j++ {
			score := self.scoreMatrix[tTypedPos{tPos{i, j}, GENERAL}]
			fmt.Printf("%4d ", int(score))
			if score >= maxScore {
				maxScore = score
				endPos = tTypedPos{tPos{i, j}, GENERAL}
			}
		}
		print("\n")
	}
	var prevPos tTypedPos
	for curPos := endPos; curPos != prevPos; prevPos, curPos = curPos, self.directionMatrix[curPos] {
		fmt.Printf("%d %d %s %f\n", curPos.n, curPos.a,
			scoreTypeToString(curPos.scoreType), self.scoreMatrix[curPos])
	}

	/*for n > -1 && a > -1 {
		score := self.scoreMatrix[tTypedPos{tPos{n, a}, GENERAL}]
		fmt.Printf("%d %d %f\n", n, a, score)
		result = append(result, tScoredPos{tPos{n, a}, score})
		prevScore1, prevScore2, prevScore3 := math.Inf(-1), math.Inf(-1), math.Inf(-1)
		if n > 1 {
			prevScore1 = self.scoreMatrix[tTypedPos{tPos{n - 1, a}, GENERAL}]
		}
		if n > 1 && a > 1 {
			prevScore2 = self.scoreMatrix[tTypedPos{tPos{n - 1, a - 1}, GENERAL}]
		}
		if a > 1 {
			prevScore3 = self.scoreMatrix[tTypedPos{tPos{n, a - 1}, GENERAL}]
		}
		fmt.Printf("%f, %f, %f\n", prevScore1, prevScore2, prevScore3)
		if prevScore1*2 >= (prevScore2 + prevScore3) {
			n--
		} else if prevScore2*2 >= (prevScore1 + prevScore3) {
			n--
			a--
		} else {
			a--
		}
	}*/
	return result
}

func (self *Alignment) getCachedScore(scoreType tScoreType, pos tPos) (float64, bool) {
	score, present := self.scoreMatrix[tTypedPos{pos, scoreType}]
	return score, present
}

func (self *Alignment) setCachedScore(scoreType tScoreType, pos tPos, score float64, prevPos tTypedPos) {
	typedPos := tTypedPos{pos, scoreType}
	self.scoreMatrix[typedPos] = score
	self.directionMatrix[typedPos] = prevPos
}

func (self *Alignment) nSeqLen() int {
	return len(self.nSeq)
}

func (self *Alignment) aSeqLen() int {
	return len(self.aSeq)
}

func (self *Alignment) getNA(nPos int) n.NucleicAcid {
	return self.nSeq[nPos-1]
}

func (self *Alignment) getAA(aPos int) a.AminoAcid {
	return self.aSeq[aPos-1]
}

func (self *Alignment) calcExtInsScore(pos tPos) float64 {
	score, present := self.getCachedScore(EXT_INS, pos)
	if present {
		return score
	}
	q := self.gapOpenPenalty
	r := self.gapExtensionPenalty
	var prevPos tTypedPos
	if pos.n == 0 && pos.a > 0 {
		score, prevPos = CalcGapScore(0, q, r), tTypedPos{tPos{0, 0}, GENERAL}
	} else {
		funcChan := make(chan func() (float64, tTypedPos))
		go func() {
			gapPenalty1 := .0
			gapPenalty2 := .0
			if pos.a > 0 && pos.a < self.aSeqLen() {
				gapPenalty1 = CalcGapScore(1, 0, r)
				gapPenalty2 = CalcGapScore(1, q, r)
			}
			prevPos := posShift(pos, -1, 0)
			funcChan <- func() (float64, tTypedPos) {
				return self.calcExtInsScore(prevPos) + gapPenalty1,
					tTypedPos{prevPos, EXT_INS}
			}
			funcChan <- func() (float64, tTypedPos) {
				return self.calcScore(prevPos) + gapPenalty2,
					tTypedPos{prevPos, GENERAL}
			}
			close(funcChan)
		}()
		score, prevPos = flowMax(funcChan)
	}
	self.setCachedScore(EXT_INS, pos, score, prevPos)
	return score
}

func (self *Alignment) calcDelScore(pos tPos) float64 {
	score, present := self.getCachedScore(DEL, pos)
	if present {
		return score
	}
	q := self.gapOpenPenalty
	r := self.gapExtensionPenalty
	var prevPos tTypedPos
	if pos.n > 0 && pos.a == 0 {
		score, prevPos = CalcGapScore(0, q, r), tTypedPos{tPos{0, 0}, GENERAL}
	} else {
		funcChan := make(chan func() (float64, tTypedPos))
		curNA := self.getNA(pos.n)
		var prevNA n.NucleicAcid
		curAA := self.getAA(pos.a)
		pos01 := posShift(pos, 0, -1)
		pos11 := posShift(pos, -1, -1)
		pos21 := posShift(pos, -2, -1)
		if pos.n < self.nSeqLen() {
			go func() {
				funcChan <- func() (float64, tTypedPos) {
					return self.calcDelScore(pos01) +
						CalcGapScore(3, 0, r), tTypedPos{pos01, DEL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos01) +
						CalcGapScore(3, q, r), tTypedPos{pos01, GENERAL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos11) +
						CalcMutationScore(curNA, n.N, n.N, curAA) +
						CalcGapScore(2, q, r), tTypedPos{pos11, GENERAL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos11) +
						CalcMutationScore(n.N, curNA, n.N, curAA) +
						CalcGapScore(1, q, r)*2, tTypedPos{pos11, GENERAL}
				}
				if pos.n > 1 {
					prevNA = self.getNA(pos.n - 1)
					funcChan <- func() (float64, tTypedPos) {
						return self.calcDelScore(pos11) +
							CalcMutationScore(n.N, curNA, n.N, curAA) +
							CalcGapScore(2, q, r), tTypedPos{pos11, DEL}
					}
					funcChan <- func() (float64, tTypedPos) {
						return self.calcScore(pos21) +
							CalcMutationScore(prevNA, curNA, n.N, curAA) +
							CalcGapScore(1, q, r), tTypedPos{pos21, GENERAL}
					}
				}
				close(funcChan)
			}()
		} else {
			go func() {
				funcChan <- func() (float64, tTypedPos) {
					return self.calcDelScore(pos01), tTypedPos{pos01, DEL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos01), tTypedPos{pos01, GENERAL}
				}
				close(funcChan)
			}()
		}
		score, prevPos = flowMax(funcChan)
	}
	self.setCachedScore(DEL, pos, score, prevPos)
	return score
}

func (self *Alignment) calcScore(pos tPos) float64 {
	score, present := self.getCachedScore(GENERAL, pos)
	if present {
		return score
	}
	var prevPos tTypedPos
	if pos.n == 0 || pos.a == 0 {
		score, prevPos = 0, tTypedPos{pos, GENERAL}
	} else {
		q := self.gapOpenPenalty
		r := self.gapExtensionPenalty
		funcChan := make(chan func() (float64, tTypedPos))
		curNA := self.getNA(pos.n)
		var prevNA, prevNA2, prevNA3 n.NucleicAcid
		curAA := self.getAA(pos.a)
		go func() {
			pos11 := posShift(pos, -1, -1)
			pos21 := posShift(pos, -2, -1)
			pos31 := posShift(pos, -3, -1)
			pos41 := posShift(pos, -4, -1)
			funcChan <- func() (float64, tTypedPos) {
				return self.calcScore(pos11) +
					CalcMutationScore(n.N, n.N, curNA, curAA) +
					CalcGapScore(2, q, r), tTypedPos{pos11, GENERAL}
			}
			if pos.n > 1 {
				prevNA = self.getNA(pos.n - 1)
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos21) +
						CalcMutationScore(prevNA, n.N, curNA, curAA) +
						CalcGapScore(1, q, r), tTypedPos{pos21, GENERAL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos21) +
						CalcMutationScore(n.N, prevNA, curNA, curAA) +
						CalcGapScore(1, q, r), tTypedPos{pos21, GENERAL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcDelScore(pos11) +
						CalcMutationScore(n.N, n.N, curNA, curAA) +
						CalcGapScore(1, 0, r), tTypedPos{pos11, DEL}
				}
			}
			if pos.n > 2 {
				prevNA2 = self.getNA(pos.n - 2)
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos31) +
							CalcMutationScore(prevNA2, prevNA, curNA, curAA),
						tTypedPos{pos31, GENERAL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcDelScore(pos21) +
						CalcMutationScore(n.N, prevNA, curNA, curAA) +
						CalcGapScore(1, 0, r), tTypedPos{pos21, DEL}
				}
			}
			if pos.n > 3 {
				prevNA3 = self.getNA(pos.n - 3)
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos41) +
						CalcMutationScore(prevNA3, prevNA, curNA, curAA) +
						CalcGapScore(1, q, r), tTypedPos{pos41, GENERAL}
				}
				funcChan <- func() (float64, tTypedPos) {
					return self.calcScore(pos41) +
						CalcMutationScore(prevNA3, prevNA2, curNA, curAA) +
						CalcGapScore(1, q, r), tTypedPos{pos41, GENERAL}
				}
			}
			funcChan <- func() (float64, tTypedPos) {
				return self.calcExtInsScore(pos), tTypedPos{pos, EXT_INS}
			}
			funcChan <- func() (float64, tTypedPos) {
				return self.calcDelScore(pos), tTypedPos{pos, DEL}
			}
			// TODO: need CODON_INS
			close(funcChan)
		}()
		score, prevPos = flowMax(funcChan)
	}
	self.setCachedScore(GENERAL, pos, score, prevPos)
	return score
}

func (self *Alignment) CalcScore() float64 {
	return self.calcScore(tPos{self.nSeqLen(), self.aSeqLen()})
}
