package cmd

import (
	"fmt"
	ap "github.com/hivdb/nucamino/alignmentprofile"
	"github.com/hivdb/nucamino/cli"
	"github.com/pkg/profile"
	"github.com/spf13/cobra"
	"io/ioutil"
	"strings"
)

// The cobra cli library will populate these variables with values
// provided as command line flags.
var alignWithInputFilename, alignWithOutputFilename, alignWithOutputFormat string
var alignWithQuiet, alignWithPprof bool
var alignWithGoroutines int

func init() {
	rootCmd.AddCommand(alignWithCmd)
	alignWithCmd.Flags().StringVarP(
		&alignWithInputFilename,
		"input-file",
		"i",
		"-",
		"input file",
	)
	alignWithCmd.Flags().StringVarP(
		&alignWithOutputFilename,
		"output-file",
		"o",
		"-",
		"output File",
	)
	alignWithCmd.Flags().StringVarP(
		&alignWithOutputFormat,
		"output-format",
		"f",
		"tsv",
		"output format. (options: \"tsv\", \"json\")",
	)
	alignWithCmd.Flags().BoolVarP(
		&alignWithQuiet,
		"quiet",
		"q",
		false,
		"hide non-error output message",
	)
	alignWithCmd.Flags().BoolVarP(
		&alignWithPprof,
		"pprof",
		"p",
		false,
		"save profiling information",
	)
	alignWithCmd.Flags().IntVar(
		&alignWithGoroutines,
		"goroutines",
		0,
		"number of goroutines the aligner will use. (default: number of CPUs)",
	)
}

func alignWithGetParameters(args []string) (*ap.AlignmentProfile, []string, error) {

	profileFileName := args[0]
	srcBytes, err := ioutil.ReadFile(profileFileName)
	if err != nil {
		return nil, nil, err
	}
	profile, err := ap.Parse(string(srcBytes))
	if err != nil {
		return nil, nil, err
	}

	genes := strings.Split(args[1], ",")
	profileGenes := profile.Genes()
	for _, gene := range genes {
		if !geneInGenes(gene, profileGenes) {
			tmpl := "%v is not an available gene in the profile %v (available genes: %v)"
			err := fmt.Errorf(tmpl, gene, profileFileName, profileGenes)
			return nil, nil, err
		}
	}

	return profile, genes, nil
}

func alignWithRun(cmd *cobra.Command, args []string) error {
	if alignWithPprof {
		defer profile.Start(profile.CPUProfile).Stop()
	}
	profile, genes, err := alignWithGetParameters(args)
	if err != nil {
		return err
	}
	return cli.PerformAlignment(
		alignWithInputFilename,
		alignWithOutputFilename,
		alignWithOutputFormat,
		genes,
		alignWithGoroutines,
		alignWithQuiet,
		*profile,
	)
}

var alignWithLongMsg = `
Loads nucleotide sequences from a FASTA file and alignts them using a
custom profile loaded from a YAML file. The first argument is the path
to the YAML file containing the profile. The second argument is a
comma separated list of genes to align against. (This list should
either be surrounded by quote marks or contain no spaces).

Examples:

	nucamino align-with custom-profile.yaml ns3
	nucamino align-with my-profile.yaml POL,GAG

You can use 'nucamino profile print' to see examples of alignment
profiles, and 'nucamino profile check' to verify that a file
represents an alignment profile that nucamino can load.

Use 'nucamino align' to use a built-in alignment profile.`

var alignWithCmd = &cobra.Command{
	Use:   "align-with",
	Short: "align sequences in a FASTA file using a custom alignment profile",
	Long:  alignWithLongMsg,
	Args:  cobra.ExactArgs(2),
	RunE:  alignWithRun,
}
