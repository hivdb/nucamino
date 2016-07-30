package scorehandler

import (
	a "../types/amino"
	n "../types/nucleic"
)

const (
	GAPINS = true
	GAPDEL = false
)

type ScoreHandler interface {
	GetSubstitutionScore(
		/* refPosition */ int,
		/* base1 */ n.NucleicAcid,
		/* base2 */ n.NucleicAcid,
		/* base3 */ n.NucleicAcid,
		/* ref */ a.AminoAcid) int
	GetGapExtensionScore( /* refPosition */ int /* isInsertion */, bool) int
	GetGapOpeningScore( /* refPosition */ int /* isINsertion */, bool) int
}
