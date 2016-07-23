package alignment

import c "../types/codon"
import a "../types/amino"
import n "../types/nucleic"
import . "../data"

type Mutation struct {
	codon     c.Codon
	consensus a.AminoAcid
}

var mutationScoreTable = make(map[Mutation]float64)

func CalcGapScore(
	gapLength int,
	gapOpenPenalty int,
	gapExtensionPenalty int) float64 {
	return float64(-gapOpenPenalty - gapLength*gapExtensionPenalty)
}

func calcMutationScore(codon c.Codon, consensus a.AminoAcid) float64 {
	// Use from cache
	mut := Mutation{codon, consensus}
	score, present := mutationScoreTable[mut]
	if present {
		return score
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
	score = float64(scores) / float64(numAAs)
	mutationScoreTable[mut] = score
	return score
}

func CalcMutationScore(
	base1 n.NucleicAcid, base2 n.NucleicAcid,
	base3 n.NucleicAcid, consensus a.AminoAcid) float64 {
	base1 = base1
	return calcMutationScore(
		c.Codon{base1, base2, base3},
		consensus)
}

func GetMutationScoreTable() map[Mutation]float64 {
	if len(mutationScoreTable) == 0 {
		for codon := range c.GenAllCodons() {
			for _, cons := range a.AminoAcids {
				mut := Mutation{codon, cons}
				mutationScoreTable[mut] = calcMutationScore(codon, cons)
			}
		}
	}
	return mutationScoreTable
}
