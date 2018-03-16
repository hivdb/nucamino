package alignmentprofile

// This package contains procedures for serializing
// rawAlignmentProfile structs. It uses the built-in template package
// instead of YAML serialization because it needs tighter control of
// the resulting data's formatting than our YAML parser has.

import (
	"bytes"
	"text/template"
)

var profileTemplateSrc = `StopCodonPenalty: {{.StopCodonPenalty}}
GapOpeningPenalty: {{.GapOpeningPenalty}}
GapExtensionPenalty: {{.GapExtensionPenalty}}
IndelCodonOpeningBonus: {{.IndelCodonOpeningBonus}}
IndelCodonExtensionBonus: {{.IndelCodonExtensionBonus}}
ReferenceSequences:
{{ range $gene, $seq := .ReferenceSequences }}  {{$gene}}:
    {{$seq}}
{{end -}}
{{ if .RawIndelScores -}} PositionalIndelScores: {{- end }}
{{range $gene, $rawIndels := .RawIndelScores}}  {{$gene}}:
{{- range $rawIndels}}
    - [ {{.Kind}}, {{.Position}}, {{.Open}}, {{.Extend}} ]
{{- end}}
{{end}}`

var profileTemplate *template.Template

func init() {
	profileTemplate = template.Must(template.New("alignmentprofile").Parse(profileTemplateSrc))
}

func Format(ap AlignmentProfile) string {
	rawProfile := ap.asRaw()
	var buff bytes.Buffer
	err := profileTemplate.Execute(&buff, rawProfile)
	if err != nil {
		panic(err)
	}
	return buff.String()
}
