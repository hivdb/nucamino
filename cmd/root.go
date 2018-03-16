package cmd

import (
	"fmt"
	"github.com/spf13/cobra"
	"os"
)

// rootCmd represents the base command when called without any subcommands
var rootCmd = &cobra.Command{
	Use:   "nucamino",
	Short: "A codon-aware viral genome aligner",
	// TODO(nknight): Add a (Long: str) more detailed description of what
	// Nucamino is, what are its capabilities and limitations, what
	// applications it's suitable/not suitable for, and where to find
	// more information.
}

// Execute adds all child commands to the root command and sets flags appropriately.
// This is called by main.main(). It only needs to happen once to the rootCmd.
func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}
