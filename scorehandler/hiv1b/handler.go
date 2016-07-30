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
	blosum62Handler *b62.Blosum62ScoreHandler
}

func (self *HIV1BScoreHandler) GetSubstitutionScore(
	position int,
	base1 n.NucleicAcid,
	base2 n.NucleicAcid,
	base3 n.NucleicAcid,
	ref a.AminoAcid) int {
	sh := self.blosum62Handler
	score := (*sh).GetSubstitutionScore(position, base1, base2, base3, ref)
	return score
}

func (self *HIV1BScoreHandler) GetGapOpeningScore(
	position int, isInsertion bool) int {
	score := self.blosum62Handler.GetGapOpeningScore(position, isInsertion)
	if self.gene == RT {
		if (isInsertion && position == 69) ||
			(!isInsertion && position >= 67 && position <= 70) {
			score = 0 // no penalty
		}
	}
	return score
}

func (self *HIV1BScoreHandler) GetGapExtensionScore(
	position int, isInsertion bool) int {
	score := self.blosum62Handler.GetGapExtensionScore(position, isInsertion)
	if self.gene == RT {
		if (isInsertion && position == 69) ||
			(!isInsertion && position >= 67 && position <= 70) {
			score = 0 // no penalty
		}
	}
	return score
}

func NewAsScoreHandler(
	gene Gene,
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int) ScoreHandler {
	return New(gene, stopCodonPenalty, gapOpenPenalty, gapExtensionPenalty, 100)
}

func New(
	gene Gene,
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int,
	scoreScale int) *HIV1BScoreHandler {
	b62handler := b62.New(stopCodonPenalty, gapOpenPenalty, gapExtensionPenalty, scoreScale)
	return &HIV1BScoreHandler{
		gene:            gene,
		scoreScale:      scoreScale,
		blosum62Handler: b62handler,
	}
}
