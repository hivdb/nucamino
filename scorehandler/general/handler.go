package general

import (
	// "fmt"
	d "github.com/hivdb/nucamino/data"
	a "github.com/hivdb/nucamino/types/amino"
	c "github.com/hivdb/nucamino/types/codon"
	n "github.com/hivdb/nucamino/types/nucleic"
)

const negInf = -int((^uint(0))>>1) - 1

type GeneralScoreHandler struct {
	scoreScale                       int
	stopCodonPenalty                 int
	gapOpenPenalty                   int
	gapExtensionPenalty              int
	indelCodonOpeningBonus           int
	indelCodonExtensionBonus         int
	positionalIndelScores            map[int]int
	positionalIndelScoresBloomFilter int
	isPositionalIndelScoreSupported  bool
	scoreMatrix                      *[a.NumAminoAcids][n.NumNucleicAcids][n.NumNucleicAcids][n.NumNucleicAcids]int
}

func (self *GeneralScoreHandler) GetCachedSubstitutionScore(
	position int,
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

func (self *GeneralScoreHandler) GetSubstitutionScoreNoCache(
	position int,
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
		// if score > 0 {
		// 	print(codon.ToString(), " ", score, "\n")
		// }
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

func (self *GeneralScoreHandler) GetSubstitutionScore(
	position int,
	base1 n.NucleicAcid,
	base2 n.NucleicAcid,
	base3 n.NucleicAcid,
	ref a.AminoAcid) int {
	// We use GetCachedSubstitutionScore separately to trigger gc's function inlining
	// which can improve the performance by about 1/6
	if score, present := self.GetCachedSubstitutionScore(position, base1, base2, base3, ref); present {
		return score
	}
	return self.GetSubstitutionScoreNoCache(position, base1, base2, base3, ref)
}

func (self *GeneralScoreHandler) GetGapOpeningScore() int {
	return -self.gapOpenPenalty
}

func (self *GeneralScoreHandler) GetGapExtensionScore() int {
	return -self.gapExtensionPenalty
}

func (self *GeneralScoreHandler) IsPositionalIndelScoreSupported() bool {
	return self.isPositionalIndelScoreSupported
}

func (self *GeneralScoreHandler) GetConstantIndelCodonScore() (int, int) {
	return self.indelCodonOpeningBonus, self.indelCodonExtensionBonus
}

func (self *GeneralScoreHandler) GetPositionalIndelCodonScore(position int, isInsertion bool) (int, int) {
	score := self.indelCodonOpeningBonus
	if self.isPositionalIndelScoreSupported {
		sign := -1 // deletion
		if isInsertion {
			sign = 1
		}
		key := sign * position
		if self.positionalIndelScoresBloomFilter&key == key {
			// the bloom filter removed most negatives; now search the real map
			_score, ok := self.positionalIndelScores[key]
			if ok {
				score = _score
			}
		}
	}
	return score, self.indelCodonExtensionBonus
}

func New(
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int,
	indelCodonOpeningBonus int,
	indelCodonExtensionBonus int,
	positionalIndelScores map[int]int,
	isPositionalIndelScoreSupported bool) *GeneralScoreHandler {
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
	positionalIndelScoresBloomFilter := 0
	for key, _ := range positionalIndelScores {
		// create a tiny bloom filter to prune negatives
		positionalIndelScoresBloomFilter |= key
	}
	return &GeneralScoreHandler{
		scoreScale:                       scoreScale,
		stopCodonPenalty:                 stopCodonPenalty * scoreScale,
		gapOpenPenalty:                   gapOpenPenalty * scoreScale,
		gapExtensionPenalty:              gapExtensionPenalty * scoreScale,
		indelCodonOpeningBonus:           indelCodonOpeningBonus * scoreScale,
		indelCodonExtensionBonus:         indelCodonExtensionBonus * scoreScale,
		positionalIndelScores:            positionalIndelScores,
		positionalIndelScoresBloomFilter: positionalIndelScoresBloomFilter,
		isPositionalIndelScoreSupported:  isPositionalIndelScoreSupported,
		scoreMatrix:                      &scoreMatrix,
	}
}
