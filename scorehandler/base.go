package scorehandler

import (
	a "../types/amino"
	n "../types/nucleic"
)

type ScoreHandler interface {
	GetSubstitutionScore(
		/* refPosition */ int,
		/* base1 */ n.NucleicAcid,
		/* base2 */ n.NucleicAcid,
		/* base3 */ n.NucleicAcid,
		/* ref */ a.AminoAcid) int
	GetInsertionScore( /* refPosition */ int) int
	GetDeletionScore( /* refPosition */ int) int
	GetGapOpeningScore( /* refPosition */ int) int
}
