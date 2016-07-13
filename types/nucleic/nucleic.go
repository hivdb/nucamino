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

var NucleicAcids = []NucleicAcid{
	A, C, G, T, W, S, M, K, R, Y, B, D, H, V, N,
}

var nucleicAcidLookup = map[NucleicAcid]rune{
	A: 'A',
	C: 'C',
	G: 'G',
	T: 'T',
	W: 'W',
	S: 'S',
	M: 'M',
	K: 'K',
	R: 'R',
	Y: 'Y',
	B: 'B',
	D: 'D',
	H: 'H',
	V: 'V',
	N: 'N',
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

func ToRune(na NucleicAcid) rune {
	return nucleicAcidLookup[na]
}

func GetUnambiguousNucleicAcids(na NucleicAcid) []NucleicAcid {
	return ambiguousNucleicAcids[na]
}
