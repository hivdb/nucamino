package alignment

import (
	. "../data"
	a "../types/amino"
	c "../types/codon"
	n "../types/nucleic"
)

const scoreScale = 100

const numN = n.NumNucleicAcids
const cubeNumN = numN * numN * numN

var mutationScoreArr [cubeNumN * a.NumAminoAcids]int

func getMutationUniq(codon c.Codon, consensus a.AminoAcid) int {
	r := (int(codon.Base1)*numN+int(codon.Base2))*numN +
		int(codon.Base3) + cubeNumN*int(consensus)
	return r
}

func CalcGapScore(
	gapLength int,
	gapOpenPenalty int,
	gapExtensionPenalty int) int {
	return (-gapOpenPenalty - gapLength*gapExtensionPenalty) * scoreScale
}

func (self *Alignment) CalcMutationScore(
	base1 n.NucleicAcid, base2 n.NucleicAcid,
	base3 n.NucleicAcid, consensus a.AminoAcid) int {
	codon := c.Codon{base1, base2, base3}
	mutUniq := getMutationUniq(codon, consensus)
	score := mutationScoreArr[mutUniq]
	if score != negInf {
		return score
	}
	scores := 0
	numAAs := 0
	aa, present := c.CodonToAminoAcid(codon)
	if present {
		score = int(LookupBlosum62(aa, consensus)) * scoreScale
	} else if isStopCodon := c.StopCodons[codon]; isStopCodon {
		score = -self.stopCodonPenalty * scoreScale
	} else {
		for unambiguousCodon := range c.GetUnambiguousCodons(codon) {
			aa, _ = c.CodonToAminoAcid(unambiguousCodon)
			scores += int(LookupBlosum62(aa, consensus))
			numAAs++
		}
		score = scores * scoreScale / numAAs
	}
	mutationScoreArr[mutUniq] = score
	return score
}

func init() {
	arrLen := cubeNumN * a.NumAminoAcids
	for i := 0; i < arrLen; i++ {
		mutationScoreArr[i] = negInf
	}
	/*for codon := range c.GenAllCodons() {
		for _, cons := range a.AminoAcids {
			calcMutationScore(codon, cons)
		}
	}*/
}
