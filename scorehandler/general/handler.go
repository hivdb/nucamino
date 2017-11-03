package general

import (
	d "github.com/hivdb/nucamino/data"
	a "github.com/hivdb/nucamino/types/amino"
	c "github.com/hivdb/nucamino/types/codon"
	n "github.com/hivdb/nucamino/types/nucleic"
	u "github.com/hivdb/nucamino/utils"
)

const posInf = int((^uint(0)) >> 1)
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
	maxPositiveScore                 int
	minPositiveScore                 int
	maxNegativeScore                 int
	scoreMatrix                      *[a.NumAminoAcids][n.NumNucleicAcids][n.NumNucleicAcids][n.NumNucleicAcids]int
}

func (self *GeneralScoreHandler) GetCachedSubstitutionScore(
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

func (self *GeneralScoreHandler) GetMaxPositiveScore() int {
	return self.maxPositiveScore
}

func (self *GeneralScoreHandler) GetMinPositiveScore() int {
	return self.minPositiveScore
}

func (self *GeneralScoreHandler) GetMaxNegativeScore() int {
	return self.maxNegativeScore
}

func (self *GeneralScoreHandler) GetSubstitutionScoreNoCache(
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
		hashed := u.SimpleFNV1a(key)
		if self.positionalIndelScoresBloomFilter&hashed == hashed {
			// the bloom filter removed most negatives; now search the real map
			_score, ok := self.positionalIndelScores[key]
			if ok {
				score = _score
			}
		}
	}
	return score, self.indelCodonExtensionBonus
}

func (self *GeneralScoreHandler) GetIndelCodonOpeningScore(position int, indelSign int) int {
	score := self.indelCodonOpeningBonus
	if self.isPositionalIndelScoreSupported {
		key := indelSign * position
		hashed := u.SimpleFNV1a(key)
		if self.positionalIndelScoresBloomFilter&hashed == hashed {
			// the bloom filter removed most negatives; now search the real map
			_score, ok := self.positionalIndelScores[key]
			if ok {
				score = _score
			}
		}
	}
	return score
}

func (self *GeneralScoreHandler) GetIndelCodonExtensionScore() int {
	return self.indelCodonExtensionBonus
}

func New(
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int,
	indelCodonOpeningBonus int,
	indelCodonExtensionBonus int,
	positionalIndelScores map[int]int,
	isPositionalIndelScoreSupported bool) *GeneralScoreHandler {
	var (
		scoreScale                         = 100
		minPositiveScore, maxNegativeScore = posInf, negInf
		maxPositiveScore                   = negInf
		scoreMatrix                        = [a.NumAminoAcids][n.NumNucleicAcids][n.NumNucleicAcids][n.NumNucleicAcids]int{}
	)

	for _, aa1 := range a.AminoAcids {
		var score int
		max, min := negInf, posInf
		for _, aa2 := range a.AminoAcids {
			score = int(d.LookupBlosum62(aa1, aa2)) * scoreScale
			if score > max {
				max = score
			}
			if score < min {
				min = score
			}
		}
		if maxPositiveScore < max {
			maxPositiveScore = max
		}
		if minPositiveScore > max {
			minPositiveScore = max
		}
		if maxNegativeScore < min {
			maxNegativeScore = min
		}
	}

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
	scaledPositionalIndelScores := map[int]int{}
	for key, score := range positionalIndelScores {
		// create a tiny bloom filter to prune negatives
		positionalIndelScoresBloomFilter |= u.SimpleFNV1a(key)
		scaledPositionalIndelScores[key] = score * scoreScale
	}
	return &GeneralScoreHandler{
		scoreScale:                       scoreScale,
		stopCodonPenalty:                 stopCodonPenalty * scoreScale,
		gapOpenPenalty:                   gapOpenPenalty * scoreScale,
		gapExtensionPenalty:              gapExtensionPenalty * scoreScale,
		indelCodonOpeningBonus:           indelCodonOpeningBonus * scoreScale,
		indelCodonExtensionBonus:         indelCodonExtensionBonus * scoreScale,
		positionalIndelScores:            scaledPositionalIndelScores,
		positionalIndelScoresBloomFilter: positionalIndelScoresBloomFilter,
		isPositionalIndelScoreSupported:  isPositionalIndelScoreSupported,
		maxPositiveScore:                 maxPositiveScore,
		minPositiveScore:                 minPositiveScore,
		maxNegativeScore:                 maxNegativeScore,
		scoreMatrix:                      &scoreMatrix,
	}
}
