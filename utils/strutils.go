package utils

import (
	"bytes"
	"unicode"
)

func StripWhiteSpace(text string) string {
	var buffer bytes.Buffer
	for _, r := range []rune(text) {
		if !unicode.IsSpace(r) {
			buffer.WriteRune(r)
		}
	}
	return buffer.String()
}
