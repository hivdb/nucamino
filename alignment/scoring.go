package alignment

import (
	. "../data"
	a "../types/amino"
	c "../types/codon"
	n "../types/nucleic"
)

type tScoreAndPresent struct {
	score   int
	present bool
}

const scoreScale = 100

var numN = n.NumNucleicAcids
var mutationScoreArr []tScoreAndPresent

func getMutationUniq(codon c.Codon, consensus a.AminoAcid) int {
	numN := n.NumNucleicAcids
	r := int(codon.Base1)
	r = r*numN + int(codon.Base2)
	r = r*numN + int(codon.Base3)
	r += numN * numN * numN * int(consensus)
	return r
}

func CalcGapScore(
	gapLength int,
	gapOpenPenalty int,
	gapExtensionPenalty int) int {
	return (-gapOpenPenalty - gapLength*gapExtensionPenalty) * scoreScale
}

func CalcMutationScore(
	base1 n.NucleicAcid, base2 n.NucleicAcid,
	base3 n.NucleicAcid, consensus a.AminoAcid) int {
	codon := c.Codon{base1, base2, base3}
	mutUniq := getMutationUniq(codon, consensus)
	sap := mutationScoreArr[mutUniq]
	if sap.present {
		return sap.score
	}
	scores := 0
	numAAs := 0
	aa, present := c.CodonToAminoAcid(codon)
	if present {
		scores += LookupBlosum62(aa, consensus)
		numAAs++
	} else {
		for unambiguousCodon := range c.GetUnambiguousCodons(codon) {
			aa, _ = c.CodonToAminoAcid(unambiguousCodon)
			scores += LookupBlosum62(aa, consensus)
			numAAs++
		}
	}
	score := scores * scoreScale / numAAs
	mutationScoreArr[mutUniq] = tScoreAndPresent{score, true}
	return score
}

func init() {
	arrLen := numN * numN * numN * a.NumAminoAcids
	mutationScoreArr = make([]tScoreAndPresent, arrLen)
	/*for codon := range c.GenAllCodons() {
		for _, cons := range a.AminoAcids {
			calcMutationScore(codon, cons)
		}
	}*/
}
