package cmd

import (
	"fmt"
	ap "github.com/hivdb/nucamino/alignmentprofile"
	"github.com/hivdb/nucamino/alignmentprofile/builtin"
	"github.com/hivdb/nucamino/cli"
	"github.com/pkg/profile"
	"github.com/spf13/cobra"
	"strings"
)

// The cobra cli library will populate these variables with values
// provided as command line flags.
var alignInputFilename, alignOutputFilename, alignOutputFormat string
var alignQuiet, alignPprof bool
var alignGoroutines int

func init() {
	rootCmd.AddCommand(alignCmd)

	alignCmd.Flags().StringVarP(
		&alignInputFilename,
		"input-file",
		"i",
		"-",
		"input file",
	)
	alignCmd.Flags().StringVarP(
		&alignOutputFilename,
		"output-file",
		"o",
		"-",
		"output File",
	)
	alignCmd.Flags().StringVarP(
		&alignOutputFormat,
		"output-format",
		"f",
		"tsv",
		"output format. (options: \"tsv\", \"json\")",
	)
	alignCmd.Flags().BoolVarP(
		&alignQuiet,
		"quiet",
		"q",
		false,
		"hide non-error output message",
	)
	alignCmd.Flags().BoolVarP(
		&alignPprof,
		"pprof",
		"p",
		false,
		"save profiling information",
	)
	alignCmd.Flags().IntVar(
		&alignGoroutines,
		"goroutines",
		0,
		"number of goroutines the aligner will use. (default: number of CPUs)",
	)
}

// Check that a gene-name is in a list of GEnes
func geneInGenes(geneArg string, genes []ap.Gene) bool {
	for _, g := range genes {
		if g.Matches(geneArg) {
			return true
		}
	}
	return false
}

func alignGetParameters(args []string) (*ap.AlignmentProfile, []string, error) {

	profileName := args[0]
	profile, found := builtin.Get(profileName)
	if !found {
		tmpl := `Unknown profile name: '%v'

See 'nucamino profile list' for a list of available profiles`
		err := fmt.Errorf(tmpl, profileName)
		return nil, nil, err
	}

	genes := strings.Split(args[1], ",")
	for idx, gene := range genes {
		genes[idx] = strings.ToUpper(gene)
	}
	profileGenes := profile.Genes()
	for _, gene := range genes {
		if !geneInGenes(gene, profileGenes) {
			tmpl := "%v is not an available gene in the profile %v (available genes: %v)"
			err := fmt.Errorf(tmpl, gene, profileName, profileGenes)
			return nil, nil, err
		}
	}

	return profile, genes, nil
}

func alignRun(cmd *cobra.Command, args []string) error {
	if alignPprof {
		defer profile.Start(profile.CPUProfile).Stop()
	}
	profile, genes, err := alignGetParameters(args)
	if err != nil {
		return err
	}
	return cli.PerformAlignment(
		alignInputFilename,
		alignOutputFilename,
		alignOutputFormat,
		genes,
		alignGoroutines,
		alignQuiet,
		*profile,
	)
}

var alignLongMsg = `
Loads nucleotide sequences from a FASTA file and aligns them using a
built-in profile. The first argument is the name of the built-in
profile to use for the alignment. The second argument is a comma
separated list of genes to align against. (This list should either be
surrounded by quote marks or contain no spaces).

Examples:

	nucamino align hiv1b pol
	nucamino align hcv1a NS3,NS5B
	nucamino align hiv1b 'gag, pol'

See 'nucamino profile list' for the available alignment profiles.

Use 'nucamino align-with' to use a custom alignment profile.`

var alignCmd = &cobra.Command{
	Use:   "align <profile name> <genes> [flags]",
	Short: "align sequences in a FASTA file using a built-in alignment profile",
	Long:  alignLongMsg,
	Args:  cobra.ExactArgs(2),
	RunE:  alignRun,
}
