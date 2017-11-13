package frameshift

import (
	"fmt"
	n "github.com/hivdb/nucamino/types/nucleic"
)

type Indel bool

const (
	INSERTION = true
	DELETION  = false
)

type FrameShift struct {
	Position         int
	NAPosition       int
	nas              []n.NucleicAcid
	NucleicAcidsText string
	indel            Indel
	IsInsertion      bool
	IsDeletion       bool
	GapLength        int
}

func New(
	position, naPosition int, nas []n.NucleicAcid,
	indel Indel, gapLength int) *FrameShift {
	return &FrameShift{
		Position:         position,
		NAPosition:       naPosition,
		nas:              nas,
		NucleicAcidsText: n.WriteString(nas),
		indel:            indel,
		IsInsertion:      indel == INSERTION,
		IsDeletion:       indel == DELETION,
		GapLength:        gapLength,
	}
}

func MakeFrameShift(position, naPosition int, allNAs []n.NucleicAcid) *FrameShift {
	//fmt.Printf("%d %s %s\n", position, n.WriteString(nas), a.WriteString(refs))
	lenAllNAs := len(allNAs)
	var frameshift *FrameShift
	if lenAllNAs%3 == 0 {
		frameshift = nil
	} else if lenAllNAs > 3 {
		frameshift = New(
			position,
			naPosition+lenAllNAs/3*3,
			allNAs[lenAllNAs/3*3:],
			INSERTION,
			lenAllNAs%3)
	} else if lenAllNAs > 0 {
		// XXX: this simplified different types of frameshift at same codon
		// position and treat them as same frameshift.
		// For example, codon "-A-" is treated as a single frameshift "del2bp"
		// but not two separated ones.
		frameshift = New(
			position,
			naPosition+lenAllNAs,
			nil,
			DELETION,
			3-lenAllNAs)
	}
	return frameshift
}

func (self *FrameShift) ToString() string {
	indel := "del"
	if self.IsInsertion {
		indel = "ins"
	}
	r := fmt.Sprintf("%d%s%dbp", self.Position, indel, self.GapLength)
	if self.IsInsertion {
		r += "_" + self.NucleicAcidsText
	}
	return r
}
