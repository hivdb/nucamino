package utils

import (
	"bytes"
	"strings"
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

func PadRightSpace(text string, length int) string {
	return text + strings.Repeat(" ", (length-len(text)))
}
