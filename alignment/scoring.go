package alignment

import (
	. "../data"
	a "../types/amino"
	c "../types/codon"
	n "../types/nucleic"
)

const scoreScale = 100

const numN = uint8(n.NumNucleicAcids)

var mutationScoreMatrix = [a.NumAminoAcids][numN][numN][numN]int{}

// TODO: extend this to support position
func (self *Alignment) CalcMutationScore(
	base1 n.NucleicAcid, base2 n.NucleicAcid,
	base3 n.NucleicAcid, consensus a.AminoAcid) int {
	codon := c.Codon{base1, base2, base3}
	score := mutationScoreMatrix[consensus][codon.Base1][codon.Base2][codon.Base3]
	if score != negInf {
		return score
	}
	if codon.IsAmbiguous() {
		// this loop also works with unambiguous codon,
		// but it slows thing down a little bit
		scores := 0
		numAAs := 0
		for _, aa := range codon.ToAminoAcids() {
			scores += int(LookupBlosum62(aa, consensus))
			numAAs++
		}
		score = scores * scoreScale / numAAs
	} else if isStopCodon := c.StopCodons[codon]; isStopCodon {
		score = -self.stopCodonPenalty * scoreScale
	} else {
		// unambiguous codon, can use unsafe function safely
		aa := codon.ToAminoAcidUnsafe()
		score = int(LookupBlosum62(aa, consensus)) * scoreScale
	}
	mutationScoreMatrix[consensus][codon.Base1][codon.Base2][codon.Base3] = score
	return score
}

func init() {
	for i, matrix3d := range mutationScoreMatrix {
		for x, matrix2d := range matrix3d {
			for y, matrix1d := range matrix2d {
				for z := range matrix1d {
					mutationScoreMatrix[i][x][y][z] = negInf
				}
			}
		}
	}
}
