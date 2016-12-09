package main

import (
	hiv1bcli "github.com/hivdb/nucamino/cli/hiv1b"
)

type HIV1BOptions struct {
	Genes                    []string    `short:"g" long:"gene" required:"yes" choice:"PR" choice:"RT" choice:"IN" description:"gene(s) the input sequences should be aligned with"`
	IndelCodonOpeningBonus   int         `long:"indel-codon-opening-bonus" value-name:"BONUS" description:"bonus score when a indel codon was opened" default:"0"`
	IndelCodonExtensionBonus int         `long:"indel-codon-extension-bonus" value-name:"BONUS" description:"bonus score when a indel codon was extended" default:"2"`
	StopCodonPenalty         int         `long:"stop-codon-penalty" value-name:"PENALTY" description:"penalty score when a stop codon was met" default:"4"`
	GapOpeningPenalty        int         `long:"gap-opening-penalty" value-name:"PENALTY" description:"penalty score when a gap was opened" default:"10"`
	GapExtensionPenalty      int         `long:"gap-extension-penalty" value-name:"PENALTY" description:"penalty score when a gap was extended" default:"2"`
	Files                    FileOptions `group:"File Options"`
}

var hiv1bOptions HIV1BOptions

func (self *HIV1BOptions) Execute(args []string) error {
	hiv1bcli.PerformAlignment(
		string(self.Files.Input), string(self.Files.Output), self.Genes,
		self.IndelCodonOpeningBonus, self.IndelCodonExtensionBonus,
		self.StopCodonPenalty, self.GapOpeningPenalty,
		self.GapExtensionPenalty)
	return nil
}

func init() {
	parser.AddCommand(
		"hiv1b", "Align HIV-1 type B sequences",
		"Use HIV-1 type B consensus from LANL to align input sequences; support genes PR, RT and IN",
		&hiv1bOptions)
}
