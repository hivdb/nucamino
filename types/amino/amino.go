package amino

import (
	"github.com/hivdb/nucamino/utils"
	"strings"
)

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

var aminoAcidLookup = [NumAminoAcids]string{
	"A",
	"C",
	"D",
	"E",
	"F",
	"G",
	"H",
	"I",
	"K",
	"L",
	"M",
	"N",
	"P",
	"Q",
	"R",
	"S",
	"T",
	"V",
	"W",
	"Y",
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

func ToString(aa AminoAcid) string {
	return aminoAcidLookup[aa]
}

func ReadString(aminoAcidSequence string) []AminoAcid {
	aminoAcidSequence = strings.ToUpper(
		utils.StripWhiteSpace(aminoAcidSequence))
	result := make([]AminoAcid, len(aminoAcidSequence))
	idx := 0
	for _, runeVal := range aminoAcidSequence {
		aa, present := aminoAcidLookupR[runeVal]
		if !present {
			continue
		}
		result[idx] = aa
		idx++
	}
	return result[:idx]
}

func WriteString(aas []AminoAcid) string {
	var result string
	for _, aa := range aas {
		result += aminoAcidLookup[aa]
	}
	return result
}
