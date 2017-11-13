package codon

import (
	a "github.com/hivdb/nucamino/types/amino"
	. "github.com/hivdb/nucamino/types/nucleic"
)

type Codon struct {
	Base1 NucleicAcid
	Base2 NucleicAcid
	Base3 NucleicAcid
}

var stopCodons = map[Codon]bool{
	Codon{T, A, A}: true,
	Codon{T, A, G}: true,
	Codon{T, G, A}: true,
	Codon{T, A, R}: true, // TAA & TAG
	Codon{T, R, A}: true, // TAA & TGA
}

var CodonToAminoAcidTable = map[Codon]a.AminoAcid{
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

var aminoAcidSearchMatrix = [a.NumAminoAcids][3][NumNucleicAcids][]Codon{}

func init() {
	for codon, aa := range CodonToAminoAcidTable {
		for idx, na := range codon.GetNucleicAcids() {
			codons := aminoAcidSearchMatrix[aa][idx][na]
			if codons == nil {
				codons = make([]Codon, 0, 4)
			}
			codons = append(codons, codon)
			aminoAcidSearchMatrix[aa][idx][na] = codons
		}
	}
}

func hasCommonCodon(codons0 []Codon, codons1 []Codon) bool {
	for _, codon0 := range codons0 {
		for _, codon1 := range codons1 {
			if codon0 == codon1 {
				return true
			}
		}
	}
	return false
}

var twoNAsCases = [3][2]int{
	[2]int{0, 1},
	[2]int{1, 2},
	[2]int{0, 2},
}

func FindBestMatch(nas []NucleicAcid, aa a.AminoAcid) Codon {
	var (
		codon, partialCodon *Codon
		lenNAs              = len(nas)
	)
	if lenNAs == 3 {
		codon = &Codon{nas[0], nas[1], nas[2]}
	} else if lenNAs == 2 {
		partialCodon = &Codon{nas[0], nas[1], N} // no match fallback
		for _, p := range twoNAsCases {
			codons0, codons1 := []Codon{}, []Codon{}
			for _, na := range GetUnambiguousNucleicAcids(nas[0]) {
				codons0 = append(codons0, aminoAcidSearchMatrix[aa][p[0]][na]...)
			}
			for _, na := range GetUnambiguousNucleicAcids(nas[1]) {
				codons1 = append(codons1, aminoAcidSearchMatrix[aa][p[1]][na]...)
			}
			pNAs := [3]NucleicAcid{N, N, N}
			pNAs[p[0]] = nas[0]
			pNAs[p[1]] = nas[1]
			pCodon := &Codon{pNAs[0], pNAs[1], pNAs[2]}
			if hasCommonCodon(codons0, codons1) {
				codon = pCodon
				break
			} else if len(codons0) > 0 || len(codons1) > 0 {
				partialCodon = pCodon
			}
		}
	} else if lenNAs == 1 {
		partialCodon = &Codon{nas[0], N, N} // no match fallback
		for i := 0; i < 3; i++ {
			codons := []Codon{}
			for _, na := range GetUnambiguousNucleicAcids(nas[0]) {
				codons = append(codons, aminoAcidSearchMatrix[aa][i][na]...)
			}
			if len(codons) > 0 {
				pNAs := [3]NucleicAcid{N, N, N}
				pNAs[i] = nas[0]
				codon = &Codon{pNAs[0], pNAs[1], pNAs[2]}
				break
			}
		}
	}
	if codon == nil {
		codon = partialCodon
	}
	return *codon
}

func (self *Codon) IsStopCodon() bool {
	return stopCodons[*self]
}

func (self *Codon) GetNucleicAcids() [3]NucleicAcid {
	return [3]NucleicAcid{
		self.Base1,
		self.Base2,
		self.Base3,
	}
}

func (self *Codon) ToAminoAcidsText() string {
	aas := make([]a.AminoAcid, 0, 1)
	hasStopCodon := false

UnambiguousCodonsLoop:
	for _, ucodon := range self.GetUnambiguousCodons() {
		if ucodon.IsStopCodon() {
			hasStopCodon = true
			continue
		}
		aa := ucodon.ToAminoAcidUnsafe()
		if len(aas) > 0 {
			for _, knownAA := range aas {
				if knownAA == aa {
					continue UnambiguousCodonsLoop
				}
			}
		}
		aas = append(aas, aa)
	}
	text := a.WriteString(aas)
	if hasStopCodon {
		text += "*"
	}
	return text
}

func (self *Codon) IsAmbiguous() bool {
	return self.Base1.IsAmbiguous() ||
		self.Base2.IsAmbiguous() ||
		self.Base3.IsAmbiguous()
}

// NOTE: This method doesn't check if self is stop codon or ambiguous
func (self *Codon) ToAminoAcidUnsafe() a.AminoAcid {
	return CodonToAminoAcidTable[*self]
}

func (self *Codon) GetUnambiguousCodons() []Codon {
	codons := make([]Codon, 0, 2)
	for _, na1 := range GetUnambiguousNucleicAcids(self.Base1) {
		for _, na2 := range GetUnambiguousNucleicAcids(self.Base2) {
			for _, na3 := range GetUnambiguousNucleicAcids(self.Base3) {
				codons = append(codons, Codon{na1, na2, na3})
			}
		}
	}
	return codons
}

func (self *Codon) ToString() string {
	return self.Base1.ToString() + self.Base2.ToString() + self.Base3.ToString()
}
