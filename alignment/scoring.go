package alignment

import c "../types/codon"
import a "../types/amino"
import . "../data"

type Mutation struct {
	codon     c.Codon
	consensus a.AminoAcid
}

var mutationScoreTable = make(map[Mutation]float64)

func GetMutationScore(codon c.Codon, consensus a.AminoAcid) float64 {
	// Use from cache
	score, present := mutationScoreTable[Mutation{codon, consensus}]
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
	return float64(scores) / float64(numAAs)
}

func GetMutationScoreTable() map[Mutation]float64 {
	if len(mutationScoreTable) == 0 {
		for codon := range c.GenAllCodons() {
			for _, cons := range a.AminoAcids {
				mut := Mutation{codon, cons}
				mutationScoreTable[mut] = GetMutationScore(codon, cons)
			}
		}
	}
	return mutationScoreTable
}
