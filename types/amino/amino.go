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

const NumAminoAcids = 20

var AminoAcids = [NumAminoAcids]AminoAcid{
	A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y,
}

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
	result := make([]AminoAcid, len(aminoAcidSequence))
	for idx, runeVal := range aminoAcidSequence {
		aa, present := aminoAcidLookupR[runeVal]
		if !present {
			continue
		}
		result[idx] = aa
	}
	return result
}

func WriteString(aas []AminoAcid) string {
	var result string
	for _, aa := range aas {
		result += aminoAcidLookup[aa]
	}
	return result
}

func ToTriplet(aa AminoAcid) string {
	return aminoAcidTripletLookup[aa]
}
