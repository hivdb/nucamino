package hiv1b

import (
	. "../../scorehandler"
	a "../../types/amino"
	n "../../types/nucleic"
	b62 "../blosum62"
)

type Gene uint8

const (
	PR Gene = iota
	RT
	IN
)

const negInf = -int((^uint(0))>>1) - 1

type HIV1BScoreHandler struct {
	gene            Gene
	scoreScale      int
	indelCodonBonus int
	blosum62Handler *b62.Blosum62ScoreHandler
}

func (self *HIV1BScoreHandler) GetSubstitutionScore(
	position int,
	base1 n.NucleicAcid,
	base2 n.NucleicAcid,
	base3 n.NucleicAcid,
	ref a.AminoAcid) int {
	bh := self.blosum62Handler
	// We use GetCachedSubstitutionScore separately to trigger gc's function inlining
	// which can improve the performance by about 1/6
	if score, present := bh.GetCachedSubstitutionScore(base1, base2, base3, ref); present {
		return score
	}
	return bh.GetSubstitutionScoreNoCache(base1, base2, base3, ref)
}

func (self *HIV1BScoreHandler) GetGapOpeningScore() int {
	return self.blosum62Handler.GetGapOpeningScore()
}

func (self *HIV1BScoreHandler) GetGapExtensionScore() int {
	return self.blosum62Handler.GetGapExtensionScore()
}

func (self *HIV1BScoreHandler) IsPositionalIndelScoreSupported() bool {
	return self.gene == RT
}

func (self *HIV1BScoreHandler) GetIndelScore(position int, isInsertion bool) int {
	score := self.indelCodonBonus
	if self.gene == RT {
		if isInsertion && position == 69 {
			score = 4 * self.scoreScale
		}
	}
	return score
}

func NewAsScoreHandler(
	gene Gene,
	indelCodonBonus int,
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int) ScoreHandler {
	return New(gene, indelCodonBonus, stopCodonPenalty, gapOpenPenalty, gapExtensionPenalty, 100)
}

func New(
	gene Gene,
	indelCodonBonus int,
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int,
	scoreScale int) *HIV1BScoreHandler {
	b62handler := b62.New(stopCodonPenalty, gapOpenPenalty, gapExtensionPenalty, scoreScale)
	return &HIV1BScoreHandler{
		gene:            gene,
		indelCodonBonus: indelCodonBonus * scoreScale,
		scoreScale:      scoreScale,
		blosum62Handler: b62handler,
	}
}
