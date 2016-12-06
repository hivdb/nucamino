package blosum62

import (
	d "github.com/hivdb/nucamino/data"
	. "github.com/hivdb/nucamino/scorehandler"
	a "github.com/hivdb/nucamino/types/amino"
	c "github.com/hivdb/nucamino/types/codon"
	n "github.com/hivdb/nucamino/types/nucleic"
)

const negInf = -int((^uint(0))>>1) - 1

type Blosum62ScoreHandler struct {
	scoreScale          int
	stopCodonPenalty    int
	gapOpenPenalty      int
	gapExtensionPenalty int
	scoreMatrix         [a.NumAminoAcids][n.NumNucleicAcids][n.NumNucleicAcids][n.NumNucleicAcids]int
}

func (self *Blosum62ScoreHandler) GetCachedSubstitutionScore(
	base1 n.NucleicAcid,
	base2 n.NucleicAcid,
	base3 n.NucleicAcid,
	ref a.AminoAcid) (int, bool) {
	score, present := self.scoreMatrix[ref][base1][base2][base3], true
	if score == negInf {
		present = false
	}
	return score, present
}

func (self *Blosum62ScoreHandler) GetSubstitutionScoreNoCache(
	base1 n.NucleicAcid,
	base2 n.NucleicAcid,
	base3 n.NucleicAcid,
	ref a.AminoAcid) (score int) {
	codon := c.Codon{base1, base2, base3}
	if codon.IsAmbiguous() {
		// this loop also works with unambiguous codon,
		// but it slows thing down a little bit
		scores := 0
		numAAs := 0
		for _, ucodon := range codon.GetUnambiguousCodons() {
			if ucodon.IsStopCodon() {
				scores -= self.stopCodonPenalty
			}
			scores += int(d.LookupBlosum62(ucodon.ToAminoAcidUnsafe(), ref)) * self.scoreScale
			numAAs++
		}
		score = scores / numAAs
		if score > 0 {
			//print(codon.ToString(), " ", score, "\n")
		}
	} else if codon.IsStopCodon() {
		score = -self.stopCodonPenalty
	} else {
		// unambiguous codon, can use unsafe function safely
		aa := codon.ToAminoAcidUnsafe()
		score = int(d.LookupBlosum62(aa, ref)) * self.scoreScale
	}
	self.scoreMatrix[ref][base1][base2][base3] = score
	return
}

func (self *Blosum62ScoreHandler) GetSubstitutionScore(
	position int,
	base1 n.NucleicAcid,
	base2 n.NucleicAcid,
	base3 n.NucleicAcid,
	ref a.AminoAcid) int {
	if score, present := self.GetCachedSubstitutionScore(base1, base2, base3, ref); present {
		return score
	}
	return self.GetSubstitutionScoreNoCache(base1, base2, base3, ref)
}

func (self *Blosum62ScoreHandler) GetGapOpeningScore() int {
	return -self.gapOpenPenalty
}

func (self *Blosum62ScoreHandler) GetGapExtensionScore() int {
	return -self.gapExtensionPenalty
}

func (self *Blosum62ScoreHandler) IsPositionalIndelScoreSupported() bool {
	return false
}

func (self *Blosum62ScoreHandler) GetConstantIndelCodonScore() (int, int) {
	return 0, 0
}

func (self *Blosum62ScoreHandler) GetPositionalIndelCodonScore(position int, isInsertion bool) (int, int) {
	return 0, 0
}

func NewAsScoreHandler(
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int) ScoreHandler {
	return New(stopCodonPenalty, gapOpenPenalty, gapExtensionPenalty, 100)
}

func New(
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int,
	scoreScale int) *Blosum62ScoreHandler {
	scoreMatrix := [a.NumAminoAcids][n.NumNucleicAcids][n.NumNucleicAcids][n.NumNucleicAcids]int{}
	for i, matrix3d := range scoreMatrix {
		for x, matrix2d := range matrix3d {
			for y, matrix1d := range matrix2d {
				for z := range matrix1d {
					scoreMatrix[i][x][y][z] = negInf
				}
			}
		}
	}
	return &Blosum62ScoreHandler{
		scoreScale:          scoreScale,
		stopCodonPenalty:    stopCodonPenalty * scoreScale,
		gapOpenPenalty:      gapOpenPenalty * scoreScale,
		gapExtensionPenalty: gapExtensionPenalty * scoreScale,
		scoreMatrix:         scoreMatrix,
	}

}
