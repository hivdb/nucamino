package cmd

import (
	"fmt"
	ap "github.com/hivdb/nucamino/alignmentprofile"
	"github.com/spf13/cobra"
	"io/ioutil"
	"os"
)

// checkCmd represents the check command
var checkCmd = &cobra.Command{
	Use:   "check [filename]",
	Short: "verify that nucamino can load the given custom alignment profile",
	Args:  cobra.RangeArgs(0, 1),
	Long:  `
Loads a YAML document and parses an alignment profile from it. Checks
that the required 'ReferenceSequences' value is present, that amino
acid sequences are valid, and that algorithm parameters have the
appropriate types. The argument, if given, is the filename to load the
profile from; reads from stdin if no argument is given.

Example:

	nucamino profile check test.yaml
	nucamino profile check < test.yaml

See 'nucamino profile print' for examples for profiles.`,
	Run:   runCheck,
}

func init() {
	profileCmd.AddCommand(checkCmd)
}

func runCheck(cmd *cobra.Command, args []string) {
	switch len(args) {
	case 0:
		checkStandardInput()
	case 1:
		checkFile(args[0])
	default:
		fmt.Fprintf(
			os.Stderr,
			"Invalid number of args: got %v, expecting 0 or 1",
			len(args),
		)
	}
}

func checkSource(source []byte) {
	_, err := ap.Parse(string(source))
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error parsing alignment profile: %v\n", err)
		os.Exit(1)
	}
	fmt.Fprintf(os.Stdout, "This profile is valid\n")
}

func checkFile(filename string) {
	srcBytes, err := ioutil.ReadFile(filename)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening file: %v\n", err)
		os.Exit(1)
	}
	checkSource(srcBytes)
}

func checkStandardInput() {
	srcBytes, err := ioutil.ReadAll(os.Stdin)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error reading from standard input: %v\n", err)
		os.Exit(1)
	}
	checkSource(srcBytes)
}
