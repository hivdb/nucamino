package nucleic

type NucleicAcid int

const (
	A NucleicAcid = iota
	C
	G
	T
	W
	S
	M
	K
	R
	Y
	B
	D
	H
	V
	N
)

const NumNucleicAcids = 15

var NucleicAcids = [NumNucleicAcids]NucleicAcid{
	A, C, G, T, W, S, M, K, R, Y, B, D, H, V, N,
}

var nucleicAcidLookup = map[NucleicAcid]string{
	A: "A",
	C: "C",
	G: "G",
	T: "T",
	W: "W",
	S: "S",
	M: "M",
	K: "K",
	R: "R",
	Y: "Y",
	B: "B",
	D: "D",
	H: "H",
	V: "V",
	N: "N",
}

var nucleicAcidLookupR = map[rune]NucleicAcid{
	'A': A,
	'C': C,
	'G': G,
	'T': T,
	'W': W,
	'S': S,
	'M': M,
	'K': K,
	'R': R,
	'Y': Y,
	'B': B,
	'D': D,
	'H': H,
	'V': V,
	'N': N,
}

var ambiguousNucleicAcids = map[NucleicAcid][]NucleicAcid{
	A: {A},
	C: {C},
	G: {G},
	T: {T},
	W: {A, T},
	S: {C, G},
	M: {A, C},
	K: {G, T},
	R: {A, G},
	Y: {C, T},
	B: {C, G, T},
	D: {A, G, T},
	H: {A, C, T},
	V: {A, C, G},
	N: {A, C, G, T},
}

func ToString(na NucleicAcid) string {
	return nucleicAcidLookup[na]
}

func ReadString(nucleicAcidSequence string) []NucleicAcid {
	result := make([]NucleicAcid, len(nucleicAcidSequence))
	for idx, runeVal := range nucleicAcidSequence {
		na, present := nucleicAcidLookupR[runeVal]
		if !present {
			na = N
		}
		result[idx] = na
	}
	return result
}

func WriteString(nas []NucleicAcid) string {
	var result string
	for _, na := range nas {
		result += nucleicAcidLookup[na]
	}
	return result
}

func GetUnambiguousNucleicAcids(na NucleicAcid) []NucleicAcid {
	return ambiguousNucleicAcids[na]
}
