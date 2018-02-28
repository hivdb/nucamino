package hcv1a

import (
	"github.com/hivdb/nucamino/alignment"
	cli "github.com/hivdb/nucamino/cli"
	data "github.com/hivdb/nucamino/data"
)

var hcv1aPositionalIndelScores = alignment.PositionalIndelScores{}

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
