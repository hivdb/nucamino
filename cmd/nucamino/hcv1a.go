package main

import (
	hcv1a "github.com/hivdb/nucamino/alignmentprofile/builtin/hcv1a"
	cli "github.com/hivdb/nucamino/cli"
	"github.com/pkg/profile"
)

type HCV1AOptions struct {
	Quiet                    bool         `short:"q" long:"quiet" description:"hide non-error information output"`
	Genes                    []string     `short:"g" long:"gene" required:"yes" choice:"NS3" choice:"NS5A" choice:"NS5B" description:"gene(s) the input sequences should be aligned with"`
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

var hcv1aOptions HCV1AOptions

func (self *HCV1AOptions) Execute(args []string) error {
	if self.Pprof.Pprof {
		defer profile.Start(profile.CPUProfile).Stop()
	}

	ioParams := cli.IOParameters{
		InputFileName:  string(self.Files.Input),
		OutputFileName: string(self.Files.Output),
		OutputFormat:   self.OutputFormat,
	}

	cli.PerformAlignment(
		ioParams.InputFileName,
		ioParams.OutputFileName,
		ioParams.OutputFormat,
		self.Genes,
		self.Goroutines,
		self.Quiet,
		hcv1a.Profile)

	return nil
}

func init() {
	parser.AddCommand(
		"hcv1a", "Align HCV-1 genotype A sequences",
		"Align input sequences to HCV-1 genotype B using the FDA recommended reference (Genebank accession no. NC_004102). Supported genes: NS3, NS5A, NS5B.",
		&hcv1aOptions)
}
