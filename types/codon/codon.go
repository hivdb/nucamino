package codon

import (
	a "../amino"
	. "../nucleic"
)

type Codon struct {
	Base1 NucleicAcid
	Base2 NucleicAcid
	Base3 NucleicAcid
}

var CodonToAminoAcidTable = map[Codon]a.AminoAcid{
	Codon{T, T, T}: a.F,
	Codon{T, T, C}: a.F,
	Codon{G, G, G}: a.G,
	Codon{T, T, T}: a.F,
	Codon{T, T, C}: a.F,
	Codon{T, T, A}: a.L,
	Codon{T, T, G}: a.L,

	Codon{C, T, T}: a.L,
	Codon{C, T, C}: a.L,
	Codon{C, T, A}: a.L,
	Codon{C, T, G}: a.L,

	Codon{A, T, T}: a.I,
	Codon{A, T, C}: a.I,
	Codon{A, T, A}: a.I,
	Codon{A, T, G}: a.M,

	Codon{G, T, T}: a.V,
	Codon{G, T, C}: a.V,
	Codon{G, T, A}: a.V,
	Codon{G, T, G}: a.V,

	Codon{T, C, T}: a.S,
	Codon{T, C, C}: a.S,
	Codon{T, C, A}: a.S,
	Codon{T, C, G}: a.S,

	Codon{C, C, T}: a.P,
	Codon{C, C, C}: a.P,
	Codon{C, C, A}: a.P,
	Codon{C, C, G}: a.P,

	Codon{A, C, T}: a.T,
	Codon{A, C, C}: a.T,
	Codon{A, C, A}: a.T,
	Codon{A, C, G}: a.T,

	Codon{G, C, T}: a.A,
	Codon{G, C, C}: a.A,
	Codon{G, C, A}: a.A,
	Codon{G, C, G}: a.A,

	Codon{T, A, T}: a.Y,
	Codon{T, A, C}: a.Y,
	//Codon{T, A, A}: a.*,
	//Codon{T, A, G}: a.*,

	Codon{C, A, T}: a.H,
	Codon{C, A, C}: a.H,
	Codon{C, A, A}: a.Q,
	Codon{C, A, G}: a.Q,

	Codon{A, A, T}: a.N,
	Codon{A, A, C}: a.N,
	Codon{A, A, A}: a.K,
	Codon{A, A, G}: a.K,

	Codon{G, A, T}: a.D,
	Codon{G, A, C}: a.D,
	Codon{G, A, A}: a.E,
	Codon{G, A, G}: a.E,

	Codon{T, G, T}: a.C,
	Codon{T, G, C}: a.C,
	//Codon{T, G, A}: a.*,//stop
	Codon{T, G, G}: a.W,

	Codon{C, G, T}: a.R,
	Codon{C, G, C}: a.R,
	Codon{C, G, A}: a.R,
	Codon{C, G, G}: a.R,

	Codon{A, G, T}: a.S,
	Codon{A, G, C}: a.S,
	Codon{A, G, A}: a.R,
	Codon{A, G, G}: a.R,

	Codon{G, G, T}: a.G,
	Codon{G, G, C}: a.G,
	Codon{G, G, A}: a.G,
	Codon{G, G, G}: a.G,
}

func CodonToAminoAcid(codon Codon) (a.AminoAcid, bool) {
	aa, present := CodonToAminoAcidTable[codon]
	return aa, present
}

func GetUnambiguousCodons(codon Codon) chan Codon {
	codonChan := make(chan Codon)
	go func() {
		for _, na1 := range GetUnambiguousNucleicAcids(codon.Base1) {
			for _, na2 := range GetUnambiguousNucleicAcids(codon.Base2) {
				for _, na3 := range GetUnambiguousNucleicAcids(codon.Base3) {
					codonChan <- Codon{na1, na2, na3}
				}
			}
		}
		close(codonChan)
	}()
	return codonChan
}

func GenAllCodons() chan Codon {
	codonChan := make(chan Codon)
	go func() {
		for _, na1 := range NucleicAcids {
			for _, na2 := range NucleicAcids {
				for _, na3 := range NucleicAcids {
					codonChan <- Codon{na1, na2, na3}
				}
			}
		}
		close(codonChan)
	}()
	return codonChan
}
