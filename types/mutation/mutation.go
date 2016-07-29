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
	insertedCodons []c.Codon
}

func New(
	position int, codon c.Codon,
	reference a.AminoAcid) *Mutation {
	return &Mutation{
		position:  position,
		codon:     &codon,
		reference: reference,
	}
}

func NewInsertion(
	position int, codon c.Codon,
	reference a.AminoAcid, insertedCodons []c.Codon) *Mutation {
	return &Mutation{
		position:       position,
		codon:          &codon,
		reference:      reference,
		isInsertion:    true,
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
	}
}

func MakeMutation(position int, nas []n.NucleicAcid, ref a.AminoAcid) *Mutation {
	lenNAs := len(nas)
	var mutation *Mutation
	if lenNAs >= 3 {
		// maybe substitution
		codon := c.Codon{nas[0], nas[1], nas[2]}
		allMatched := true
		for _, ucodon := range c.GetUnambiguousCodons(codon) {
			allMatched = allMatched && c.CodonToAminoAcidTable[ucodon] == ref
		}
		lenInsCodons := lenNAs/3 - 1
		if lenInsCodons > 0 {
			// insertion
			insertedCodons := make([]c.Codon, lenInsCodons)
			for idx := range insertedCodons {
				p := idx*3 + 3
				insertedCodons[idx] = c.Codon{nas[p], nas[p+1], nas[p+2]}
			}
			mutation = NewInsertion(position, codon, ref, insertedCodons)
		} else if !allMatched {
			mutation = New(position, codon, ref)
		}
	} else if lenNAs > 0 {
		// codon missed 1 or 2 NAs
		codon := c.FindBestMatch(nas, ref)
		mutation = New(position, codon, ref)
	} else if lenNAs == 0 {
		// deletion
		mutation = NewDeletion(position, ref)
	}
	return mutation
}

func (self *Mutation) GetPosition() int {
	return self.position
}

func (self *Mutation) GetCodon() c.Codon {
	return *self.codon
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
		r += a.WriteString(self.codon.ToAminoAcids())
		if self.isInsertion {
			r += "_"
			nas += "_"
			for _, insCodon := range self.insertedCodons {
				insertedAAs := a.WriteString(insCodon.ToAminoAcids())
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
