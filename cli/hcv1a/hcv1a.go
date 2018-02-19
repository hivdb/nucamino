package hcv1a

import (
	cli "github.com/hivdb/nucamino/cli/cli"
	data "github.com/hivdb/nucamino/data"
)

var hcv1aPositionalIndelScores = cli.PositionalIndelScores{}

func PerformAlignment(
	ioParams cli.IOParameters,
	textGenes []string,
	goroutines int,
	quiet bool,
	alignmentParams cli.AlignmentParameters) {

	cli.PerformAlignment(
		ioParams.InputFileName,
		ioParams.OutputFileName,
		ioParams.OutputFormat,
		textGenes,
		goroutines,
		quiet,
		alignmentParams,
		hcv1aPositionalIndelScores,
		data.HCV1ARefLookup)
}
