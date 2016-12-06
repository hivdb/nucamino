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
	position  int
	nas       []n.NucleicAcid
	indel     Indel
	gapLength int
}

func New(
	position int, nas []n.NucleicAcid,
	indel Indel, gapLength int) *FrameShift {
	return &FrameShift{
		position:  position,
		nas:       nas,
		indel:     indel,
		gapLength: gapLength,
	}
}

func MakeFrameShift(position int, allNAs []n.NucleicAcid) *FrameShift {
	//fmt.Printf("%d %s %s\n", position, n.WriteString(nas), a.WriteString(refs))
	lenAllNAs := len(allNAs)
	var frameshift *FrameShift
	if lenAllNAs%3 == 0 {
		frameshift = nil
	} else if lenAllNAs > 3 {
		frameshift = New(
			position,
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
			nil,
			DELETION,
			3-lenAllNAs)
	}
	return frameshift
}

func (self *FrameShift) GetPosition() int {
	return self.position
}

func (self *FrameShift) GetNucleicAcids() []n.NucleicAcid {
	return self.nas
}

func (self *FrameShift) IsInsertion() bool {
	return bool(self.indel)
}

func (self *FrameShift) IsDeletion() bool {
	return !bool(self.indel)
}

func (self *FrameShift) GetGapLength() int {
	return self.gapLength
}

func (self *FrameShift) ToString() string {
	indel := "del"
	if self.indel == INSERTION {
		indel = "ins"
	}
	r := fmt.Sprintf("%d%s%dbp", self.position, indel, self.gapLength)
	if self.indel == INSERTION {
		r += "_" + n.WriteString(self.nas)
	}
	return r
}
