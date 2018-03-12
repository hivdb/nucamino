package alignmentprofile

import (
	yaml "gopkg.in/yaml.v2"
)

func (p *AlignmentProfile) UnmarshalYAML(unmarshal func(interface{}) error) error {
	var raw rawAlignmentProfile
	err := unmarshal(&raw)
	if err != nil {
		return err
	}
	profile, err := raw.asProfile()
	if err != nil {
		return err
	}
	*p = *profile
	return nil
}

// Parse an AlignmentProfile from YAML
func Parse(src string) (*AlignmentProfile, error) {
	var profile AlignmentProfile
	err := yaml.Unmarshal([]byte(src), &profile)
	if err != nil {
		return nil, err
	}
	err = profile.validate()
	if err != nil {
		return nil, err
	}
	return &profile, nil
}
