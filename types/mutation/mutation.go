package mutation

import (
	a "../amino"
	c "../codon"
	n "../nucleic"
	"fmt"
)

type Mutation struct {
	position       int
	codon          *c.Codon
	reference      a.AminoAcid
	isInsertion    bool
	isDeletion     bool
	control        string
	insertedCodons []c.Codon
}

func New(
	position int, codon c.Codon,
	reference a.AminoAcid, control string) *Mutation {
	return &Mutation{
		position:  position,
		codon:     &codon,
		reference: reference,
		control:   control,
	}
}

func NewInsertion(
	position int, codon c.Codon, reference a.AminoAcid,
	insertedCodons []c.Codon, control string) *Mutation {
	return &Mutation{
		position:       position,
		codon:          &codon,
		reference:      reference,
		isInsertion:    true,
		control:        control,
		insertedCodons: insertedCodons,
	}
}

func NewDeletion(
	position int, reference a.AminoAcid) *Mutation {
	return &Mutation{
		position:   position,
		codon:      nil,
		reference:  reference,
		isDeletion: true,
		control:    "---",
	}
}

func MakeMutation(position int, nas []n.NucleicAcid, ref a.AminoAcid) *Mutation {
	lenNAs := len(nas)
	var (
		control  string
		mutation *Mutation
	)
	if lenNAs >= 3 {
		// maybe substitution
		codon := c.Codon{nas[0], nas[1], nas[2]}
		allMatched := true
		for _, ucodon := range codon.GetUnambiguousCodons() {
			if ucodon.IsStopCodon() {
				allMatched = false
			} else {
				allMatched = allMatched && ucodon.ToAminoAcidUnsafe() == ref
			}
		}
		if allMatched {
			control = ":::"
		} else {
			control = "..."
		}
		lenInsCodons := lenNAs/3 - 1
		if lenInsCodons > 0 {
			// insertion
			insertedCodons := make([]c.Codon, lenInsCodons)
			for idx := range insertedCodons {
				p := idx*3 + 3
				insertedCodons[idx] = c.Codon{nas[p], nas[p+1], nas[p+2]}
				control += "+++"
			}
			mutation = NewInsertion(position, codon, ref, insertedCodons, control)
		} else if !allMatched {
			mutation = New(position, codon, ref, control)
		}
	} else if lenNAs > 0 {
		// codon missed 1 or 2 NAs
		codon := c.FindBestMatch(nas, ref)
		for _, na := range codon.GetNucleicAcids() {
			if na == n.N {
				control += "-"
			} else {
				control += "."
			}
		}
		mutation = New(position, codon, ref, control)
	} else if lenNAs == 0 {
		// deletion
		mutation = NewDeletion(position, ref)
	}
	return mutation
}

func (self *Mutation) GetPosition() int {
	return self.position
}

func (self *Mutation) IsInsertion() bool {
	return self.isInsertion
}

func (self *Mutation) IsDeletion() bool {
	return self.isDeletion
}

func (self *Mutation) GetCodon() c.Codon {
	return *self.codon
}

func (self *Mutation) GetControl() string {
	return self.control
}

func (self *Mutation) GetInsertedCodons() []c.Codon {
	return self.insertedCodons
}

func (self *Mutation) GetReference() a.AminoAcid {
	return self.reference
}

func (self *Mutation) ToString() string {
	r := fmt.Sprintf("%s%d", a.ToString(self.reference), self.position)
	nas := ""
	if self.isDeletion {
		r += "-"
	} else {
		nas += ":" + self.codon.ToString()
		r += self.codon.ToAminoAcidsText()
		if self.isInsertion {
			r += "_"
			nas += "_"
			for _, insCodon := range self.insertedCodons {
				insertedAAs := insCodon.ToAminoAcidsText()
				if len(insertedAAs) > 1 {
					r += "[" + insertedAAs + "]"
				} else {
					r += insertedAAs
				}
				nas += insCodon.ToString()
			}
		}
	}
	return r + nas
}
