package blosum62

import (
	d "../../data"
	. "../../scorehandler"
	a "../../types/amino"
	c "../../types/codon"
	n "../../types/nucleic"
)

const negInf = -int((^uint(0))>>1) - 1

type tBlosum62ScoreHandler struct {
	scoreScale          int
	stopCodonPenalty    int
	gapOpenPenalty      int
	gapExtensionPenalty int
	scoreMatrix         [a.NumAminoAcids][n.NumNucleicAcids][n.NumNucleicAcids][n.NumNucleicAcids]int
}

func (self *tBlosum62ScoreHandler) GetSubstitutionScore(
	position int,
	base1 n.NucleicAcid,
	base2 n.NucleicAcid,
	base3 n.NucleicAcid,
	ref a.AminoAcid) int {
	score := self.scoreMatrix[ref][base1][base2][base3]
	if score != negInf {
		return score
	}
	codon := c.Codon{base1, base2, base3}
	if codon.IsAmbiguous() {
		// this loop also works with unambiguous codon,
		// but it slows thing down a little bit
		scores := 0
		numAAs := 0
		for _, aa := range codon.ToAminoAcids() {
			scores += int(d.LookupBlosum62(aa, ref))
			numAAs++
		}
		score = scores * self.scoreScale / numAAs
	} else if isStopCodon := c.StopCodons[codon]; isStopCodon {
		score = -self.stopCodonPenalty
	} else {
		// unambiguous codon, can use unsafe function safely
		aa := codon.ToAminoAcidUnsafe()
		score = int(d.LookupBlosum62(aa, ref)) * self.scoreScale
	}
	self.scoreMatrix[ref][base1][base2][base3] = score
	return score
}

func (self *tBlosum62ScoreHandler) GetGapOpeningScore(
	position int) int {
	return -self.gapOpenPenalty
}

func (self *tBlosum62ScoreHandler) GetInsertionScore(
	position int) int {
	return -self.gapExtensionPenalty
}

func (self *tBlosum62ScoreHandler) GetDeletionScore(
	position int) int {
	return -self.gapExtensionPenalty
}

func NewBlosum62ScoreHandler(
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int) ScoreHandler {
	scoreScale := 100
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
	return &tBlosum62ScoreHandler{
		scoreScale:          scoreScale,
		stopCodonPenalty:    stopCodonPenalty * scoreScale,
		gapOpenPenalty:      gapOpenPenalty * scoreScale,
		gapExtensionPenalty: gapExtensionPenalty * scoreScale,
		scoreMatrix:         scoreMatrix,
	}

}
