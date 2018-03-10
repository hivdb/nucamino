package cmd

import (
	"fmt"
	ap "github.com/hivdb/nucamino/alignmentprofile"
	builtin "github.com/hivdb/nucamino/alignmentprofile/builtin"
	"github.com/spf13/cobra"
	"os"
)

// printCmd represents the print command
var printCmd = &cobra.Command{
	Use:   "print <profile name>",
	Short: "Print the contents of a built-in profile",
	Args:  cobra.ExactArgs(1),
	Run: func(cmd *cobra.Command, args []string) {
		profileName := args[0]
		profile, found := builtin.Get(profileName)
		if !found {
			fmt.Fprintf(os.Stderr, "\nNo such profile built-in: %v\n\n", profileName)
			fmt.Fprintf(
				os.Stderr,
				"See the 'profile list' command for a list of available profiles.\n\n",
			)
			cmd.Usage()
			os.Exit(1)
			return
		}
		profileString := ap.Format(*profile)
		fmt.Println(profileString)
	},
}

func init() {
	profileCmd.AddCommand(printCmd)
}
