package cmd

import (
	"fmt"
	ap "github.com/hivdb/nucamino/alignmentprofile"
	"github.com/hivdb/nucamino/alignmentprofile/builtin"
	"github.com/spf13/cobra"
	"log"
	"regexp"
)

func listGenes(cmd *cobra.Command, args []string) error {
	profileName := args[0]
	profile, found := builtin.Get(profileName)
	if !found {
		tmpl := `Unknown profile name: '%v'

see 'nucamino profile list' for a list of available alignment profiles`
		err := fmt.Errorf(tmpl, profileName)
		return err
	}
	genes := profile.Genes()
	if len(args) == 2 {
		patternSrc := args[1]
		pattern, err := regexp.Compile(patternSrc)
		if err != nil {
			log.Printf("Error in search pattern")
			return err
		}
		matchingGenes := make([]ap.Gene, 0, len(genes))
		for _, gene := range genes {
			if pattern.MatchString(string(gene)) {
				matchingGenes = append(matchingGenes, gene)
			}
		}
		genes = matchingGenes
	}
	for _, gene := range(genes) {
		fmt.Println(gene)
	}
	return nil
}


var listGenesCmd = &cobra.Command{
	Use:   "list-genes profile [pattern]",
	Short: "List the available genes in a built-in alignment profile",
	Long: `This command lists the genes available  in a built-in alignment
profile. These names could be used to construct an align command, or just
to learn about the available options without printing out the whole profile.

The pattern argument is used to filter the list. It's interpreted as a
regular expression. For example:

	nucamino profile list-genes hcv1a		List the genes in the built-in HCV1a profile.
	nucamino profile list-genes	hiv1b ^G	List the genes in the HIV1b profile that start
											with a 'G'.

See https://github.com/google/re2/wiki/Syntax for the exact details of
the regular expression syntax.`,
	Args: cobra.RangeArgs(1, 2),
	RunE: listGenes,
}


func init() {
	profileCmd.AddCommand(listGenesCmd)
}
