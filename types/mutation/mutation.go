package mutation

import (
	"fmt"
	a "github.com/hivdb/nucamino/types/amino"
	c "github.com/hivdb/nucamino/types/codon"
	n "github.com/hivdb/nucamino/types/nucleic"
)

type Mutation struct {
	Position               int
	CodonText              string
	AminoAcidText          string
	codon                  *c.Codon
	ReferenceText          string
	reference              a.AminoAcid
	IsInsertion            bool
	IsDeletion             bool
	IsPartial              bool
	Control                string
	InsertedCodonsText     string
	InsertedAminoAcidsText string
	insertedCodons         []c.Codon
}

func New(
	position int, codon c.Codon,
	reference a.AminoAcid, isPartial bool, control string) *Mutation {
	codonText := []rune(codon.ToString())
	if isPartial {
		for pos, char := range control {
			if char == '-' {
				codonText[pos] = ' '
			}
		}
	}
	return &Mutation{
		Position:      position,
		codon:         &codon,
		CodonText:     string(codonText),
		AminoAcidText: codon.ToAminoAcidsText(),
		reference:     reference,
		ReferenceText: a.ToString(reference),
		IsPartial:     isPartial,
		Control:       control,
	}
}

func NewInsertion(
	position int, codon c.Codon, reference a.AminoAcid,
	insertedCodons []c.Codon, control string) *Mutation {

	var (
		insertedCodonsText     string
		insertedAminoAcidsText string
	)
	for _, insCodon := range insertedCodons {
		insertedCodonsText += insCodon.ToString()
		insertedAAs := insCodon.ToAminoAcidsText()
		if len(insertedAAs) > 1 {
			insertedAminoAcidsText += "[" + insertedAAs + "]"
		} else {
			insertedAminoAcidsText += insertedAAs
		}
	}

	mutation := New(position, codon, reference, false, control)
	mutation.IsInsertion = true
	mutation.insertedCodons = insertedCodons
	mutation.InsertedCodonsText = insertedCodonsText
	mutation.InsertedAminoAcidsText = insertedAminoAcidsText

	return mutation
}

func NewDeletion(
	position int, reference a.AminoAcid) *Mutation {
	return &Mutation{
		Position:      position,
		codon:         nil,
		CodonText:     "",
		AminoAcidText: "",
		reference:     reference,
		ReferenceText: a.ToString(reference),
		IsDeletion:    true,
		Control:       "---",
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
			mutation = New(position, codon, ref, false, control)
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
		mutation = New(position, codon, ref, true, control)
	} else if lenNAs == 0 {
		// deletion
		mutation = NewDeletion(position, ref)
	}
	return mutation
}

func (self *Mutation) GetCodon() c.Codon            { return *self.codon }
func (self *Mutation) GetReference() a.AminoAcid    { return self.reference }
func (self *Mutation) GetInsertedCodons() []c.Codon { return self.insertedCodons }

func (self *Mutation) ToString() string {
	r := fmt.Sprintf("%s%d", self.ReferenceText, self.Position)
	nas := ""
	if self.IsDeletion {
		r += "-"
	} else {
		nas += ":" + self.CodonText
		if self.IsPartial {
			r += "X" // mutation contains del gap doesn't get displayed
		} else {
			r += self.AminoAcidText
		}
		if self.IsInsertion {
			r += "_" + self.InsertedAminoAcidsText
			nas += "_" + self.InsertedCodonsText
		}
	}
	return r + nas
}
