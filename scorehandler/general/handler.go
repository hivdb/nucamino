package general

import (
	d "github.com/hivdb/nucamino/data"
	a "github.com/hivdb/nucamino/types/amino"
	c "github.com/hivdb/nucamino/types/codon"
	n "github.com/hivdb/nucamino/types/nucleic"
)

const negInf = -int((^uint(0))>>1) - 1

// FNV algorithm was used to generated index number for the bloom filter
const FNVPrime int64 = 16777619
const FNVOffsetBasis int64 = 2166136261

func simpleFNV1a(key int) int64 {
	// http://www.isthe.com/chongo/tech/comp/fnv/index.html
	// The possible values of "key" are between -50000 and 50000
	// and there's no collision in this region

	// This is not a strict implementation of FNV (the key is not an octet),
	// however it reduced the collisions of bloom filter in HIV-1 POL case
	return (FNVOffsetBasis ^ int64(key)) * FNVPrime
}

type GeneralScoreHandler struct {
	scoreScale                       int
	stopCodonPenalty                 int
	gapOpenPenalty                   int
	gapExtensionPenalty              int
	indelCodonOpeningBonus           int
	indelCodonExtensionBonus         int
	positionalIndelScores            map[int][2]int
	positionalIndelScoresBloomFilter int64
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
	ext := self.indelCodonExtensionBonus
	if self.isPositionalIndelScoreSupported {
		sign := -1 // deletion
		if isInsertion {
			sign = 1
		}
		key := sign * position
		hashed := simpleFNV1a(key)
		if self.positionalIndelScoresBloomFilter&hashed == hashed {
			// the bloom filter removed most negatives; now search the real map
			_score, ok := self.positionalIndelScores[key]
			if ok {
				score = _score[0]
				ext = _score[1]
			}
		}
	}
	return score, ext
}

func New(
	stopCodonPenalty int,
	gapOpenPenalty int,
	gapExtensionPenalty int,
	indelCodonOpeningBonus int,
	indelCodonExtensionBonus int,
	positionalIndelScores map[int][2]int,
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
	var positionalIndelScoresBloomFilter int64 = 0
	scaledPositionalIndelScores := map[int][2]int{}
	for key, score := range positionalIndelScores {
		// create a tiny bloom filter to prune negatives
		positionalIndelScoresBloomFilter |= simpleFNV1a(key)
		scaledPositionalIndelScores[key] = [2]int{score[0] * scoreScale, score[1] * scoreScale}
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
		scoreMatrix:                      &scoreMatrix,
	}
}
