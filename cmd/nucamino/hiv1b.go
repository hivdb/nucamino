package main

import (
	cli "github.com/hivdb/nucamino/cli"
	hiv1bcli "github.com/hivdb/nucamino/cli/hiv1b"
	"github.com/pkg/profile"
)

type HIV1BOptions struct {
	Quiet                    bool         `short:"q" long:"quiet" description:"hide non-error information output"`
	Genes                    []string     `short:"g" long:"gene" required:"yes" choice:"GAG" choice:"POL" choice:"GP41" description:"gene(s) the input sequences should be aligned with"`
	IndelCodonOpeningBonus   int          `long:"indel-codon-opening-bonus" value-name:"BONUS" description:"bonus score when a indel codon was opened" default:"0"`
	IndelCodonExtensionBonus int          `long:"indel-codon-extension-bonus" value-name:"BONUS" description:"bonus score when a indel codon was extended" default:"2"`
	StopCodonPenalty         int          `long:"stop-codon-penalty" value-name:"PENALTY" description:"penalty score when a stop codon was met" default:"4"`
	GapOpeningPenalty        int          `long:"gap-opening-penalty" value-name:"PENALTY" description:"penalty score when a gap was opened" default:"10"`
	GapExtensionPenalty      int          `long:"gap-extension-penalty" value-name:"PENALTY" description:"penalty score when a gap was extended" default:"2"`
	Goroutines               int          `long:"goroutines" value-name:"GOROUTINES" description:"number of goroutines the alignment will use. Use the core number when equals to 0" default:"0"`
	OutputFormat             string       `long:"output-format" value-name:"OUTPUT_FORMAT" choice:"tsv" choice:"json" description:"output format of the alignment result" default:"tsv"`
	Files                    FileOptions  `group:"File Options"`
	Pprof                    PprofOptions `group:"Pprof Options"`
}

var hiv1bOptions HIV1BOptions

func (self *HIV1BOptions) Execute(args []string) error {
	if self.Pprof.Pprof {
		defer profile.Start(profile.CPUProfile).Stop()
	}

	ioParams := cli.IOParameters{
		InputFileName:  string(self.Files.Input),
		OutputFileName: string(self.Files.Output),
		OutputFormat:   self.OutputFormat,
	}

	alignmentParams := cli.AlignmentParameters{
		IndelCodonOpeningBonus:   self.IndelCodonOpeningBonus,
		IndelCodonExtensionBonus: self.IndelCodonExtensionBonus,
		StopCodonPenalty:         self.StopCodonPenalty,
		GapOpeningPenalty:        self.GapOpeningPenalty,
		GapExtensionPenalty:      self.GapExtensionPenalty,
	}

	hiv1bcli.PerformAlignment(ioParams, self.Genes, self.Goroutines,
		self.Quiet, alignmentParams)

	return nil
}

func init() {
	parser.AddCommand(
		"hiv1b", "Align HIV-1 type B sequences",
		"Use HIV-1 type B consensus from LANL to align input sequences; support genes POL (56gag + 99PR + 560RT + 288IN)",
		&hiv1bOptions)
}
