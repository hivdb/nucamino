package frameshift

import (
	n "../nucleic"
	"fmt"
)

type Indel uint8

const (
	INSERTION Indel = iota
	DELETION
)

type FrameShift struct {
	position int
	nas      []n.NucleicAcid
	indel    Indel
	lenNAs   int
}

func New(
	position int, nas []n.NucleicAcid,
	indel Indel, lenNAs int) *FrameShift {
	return &FrameShift{
		position: position,
		nas:      nas,
		indel:    indel,
		lenNAs:   lenNAs,
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

func (self *FrameShift) ToString() string {
	indel := "del"
	if self.indel == INSERTION {
		indel = "ins"
	}
	r := fmt.Sprintf("%d%s%dbp", self.position, indel, self.lenNAs)
	if self.indel == INSERTION {
		r += "_" + n.WriteString(self.nas)
	}
	return r
}
