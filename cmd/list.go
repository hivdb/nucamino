package cmd

import (
	"fmt"
	builtin "github.com/hivdb/nucamino/alignmentprofile/builtin"
	"github.com/spf13/cobra"
	"log"
	"os"
	"regexp"
)

// listCmd represents the list command
var listCmd = &cobra.Command{
	Use:   "list [pattern]",
	Short: "List available built-in alignment profiles",
	Long: `This command lists the available built-in alignment profiles. These
names can be passed to the align command to use the built-in profiles
when aligning sequences.

The pattern argument is used to filter the list. It's interpreted as a
regular expression. For example:

    nucamino profile list '^hiv'            List all profiles that start with 'hiv'
    nucamino profile list hiv               List all profiles that contain 'hiv' anywhere
    nucamino profile list '.*1[abc]$'       List all profiles end in 1a, 1b, or 1c
    nucamino profile list 'hcv1.'           List all profiles that match 'hcv1' followed by 
                                            any other letter

See https://github.com/google/re2/wiki/Syntax for the exact details of
the regular expression syntax.`,
	Args: cobra.RangeArgs(0, 1),
	Run: func(cmd *cobra.Command, args []string) {
		profileNames := builtin.List()
		if len(args) == 1 {
			var pattern *regexp.Regexp
			patternSrc := args[0]
			pattern, err := regexp.Compile(patternSrc)
			if err != nil {
				log.Printf("Error in pattern: %v", err)
				os.Exit(1)
				return 
			}
			matchingProfileNames := make([]string, 0, len(profileNames))
			for _, name := range(profileNames) {
				if pattern.MatchString(name) {
					matchingProfileNames = append(matchingProfileNames, name)
				}
			}
			profileNames = matchingProfileNames
		}
		for _, name := range profileNames {
			fmt.Println(name)
		}
	},
}

func init() {
	profileCmd.AddCommand(listCmd)
}
