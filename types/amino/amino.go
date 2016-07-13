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

var aminoAcidLookup = map[AminoAcid]rune{
	A: 'A',
	C: 'C',
	D: 'D',
	E: 'E',
	F: 'F',
	G: 'G',
	H: 'H',
	I: 'I',
	K: 'K',
	L: 'L',
	M: 'M',
	N: 'N',
	P: 'P',
	Q: 'Q',
	R: 'R',
	S: 'S',
	T: 'T',
	V: 'V',
	W: 'W',
	Y: 'Y',
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

func ToRune(aa AminoAcid) rune {
	return aminoAcidLookup[aa]
}

func ToTriplet(aa AminoAcid) string {
	return aminoAcidTripletLookup[aa]
}
