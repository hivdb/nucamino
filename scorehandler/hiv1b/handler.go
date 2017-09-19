package hiv1b

import (
	. "github.com/hivdb/nucamino/scorehandler"
	b62 "github.com/hivdb/nucamino/scorehandler/blosum62"
	a "github.com/hivdb/nucamino/types/amino"
	n "github.com/hivdb/nucamino/types/nucleic"
)

type Gene uint8

const (
	POL Gene = iota
)

var GeneLookup = map[string]Gene{
	"POL": POL,
}

const negInf = -int((^uint(0))>>1) - 1

type HIV1BScoreHandler struct {
	gene                     Gene
	scoreScale               int
	indelCodonOpeningBonus   int
	indelCodonExtensionBonus int
	blosum62Handler          *b62.Blosum62ScoreHandler
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
	return self.gene == POL
}

func (self *HIV1BScoreHandler) GetConstantIndelCodonScore() (int, int) {
	return self.indelCodonOpeningBonus, self.indelCodonExtensionBonus
}

func (self *HIV1BScoreHandler) GetPositionalIndelCodonScore(position int, isInsertion bool) (int, int) {
	score := self.indelCodonOpeningBonus
	if self.gene == POL {
		// offset = 56prePR + 99PR = 155
		if isInsertion {
			switch position {
			case 155 + 65, 155 + 66, 155 + 67:
				score = -9 * self.scoreScale
				break
			case 155 + 63, 155 + 64, 155 + 68, 155 + 70, 155 + 71, 155 + 72, 155 + 73:
				score = -3 * self.scoreScale
				break
			case 155 + 69:
				score = 6 * self.scoreScale
				break
			}
		}
	}
	return score, self.indelCodonExtensionBonus
}

func NewAsScoreHandler(
	gene Gene,
	indelCodonOpeningBonus int,
	indelCodonExtensionBonus int,
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int) ScoreHandler {
	return New(
		gene, indelCodonOpeningBonus, indelCodonExtensionBonus,
		stopCodonPenalty, gapOpenPenalty, gapExtensionPenalty, 100)
}

func New(
	gene Gene,
	indelCodonOpeningBonus int,
	indelCodonExtensionBonus int,
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int,
	scoreScale int) *HIV1BScoreHandler {
	b62handler := b62.New(stopCodonPenalty, gapOpenPenalty, gapExtensionPenalty, scoreScale)
	return &HIV1BScoreHandler{
		gene: gene,
		indelCodonOpeningBonus:   indelCodonOpeningBonus * scoreScale,
		indelCodonExtensionBonus: indelCodonExtensionBonus * scoreScale,
		scoreScale:               scoreScale,
		blosum62Handler:          b62handler,
	}
}
