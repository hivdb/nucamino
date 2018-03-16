package cmd

import (
	"github.com/spf13/cobra"
)

// profileCmd represents the profile command
var profileCmd = &cobra.Command{
	Use:   "profile",
	Short: "Manage alignment profiles",
	Long: `An alignment profile contains all the information that nucamino needs to align
a nucleotide sequence: alignment tuning parameters, reference sequences, and the
list of custom scores for known indels.`,
}

func init() {
	rootCmd.AddCommand(profileCmd)
}
