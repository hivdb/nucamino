package amino

type AminoAcid int

const (
	A AminoAcid = iota
	C
	D
	E
	F
	G
	H
	I
	K
	L
	M
	N
	P
	Q
	R
	S
	T
	V
	W
	Y
)

var AminoAcids = []AminoAcid{
	A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
}

var NumAminoAcids = len(AminoAcids)

var aminoAcidLookup = map[AminoAcid]string{
	A: "A",
	C: "C",
	D: "D",
	E: "E",
	F: "F",
	G: "G",
	H: "H",
	I: "I",
	K: "K",
	L: "L",
	M: "M",
	N: "N",
	P: "P",
	Q: "Q",
	R: "R",
	S: "S",
	T: "T",
	V: "V",
	W: "W",
	Y: "Y",
}

var aminoAcidLookupR = map[rune]AminoAcid{
	'A': A,
	'C': C,
	'D': D,
	'E': E,
	'F': F,
	'G': G,
	'H': H,
	'I': I,
	'K': K,
	'L': L,
	'M': M,
	'N': N,
	'P': P,
	'Q': Q,
	'R': R,
	'S': S,
	'T': T,
	'V': V,
	'W': W,
	'Y': Y,
}

var aminoAcidTripletLookup = map[AminoAcid]string{
	A: "Ala",
	C: "Cys",
	D: "Asp",
	E: "Glu",
	F: "Phe",
	G: "Gly",
	H: "His",
	I: "Ile",
	K: "Lys",
	L: "Leu",
	M: "Met",
	N: "Asn",
	P: "Pro",
	Q: "Gln",
	R: "Arg",
	S: "Ser",
	T: "Thr",
	V: "Val",
	W: "Trp",
	Y: "Tyr",
}

func ToString(aa AminoAcid) string {
	return aminoAcidLookup[aa]
}

func ReadString(aminoAcidSequence string) []AminoAcid {
	result := make([]AminoAcid, 0)
	for _, runeVal := range aminoAcidSequence {
		result = append(result, aminoAcidLookupR[runeVal])
	}
	return result
}

func ToTriplet(aa AminoAcid) string {
	return aminoAcidTripletLookup[aa]
}
