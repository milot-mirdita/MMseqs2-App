package main

import (
	"crypto/sha256"
	"encoding/base64"
	"errors"
	"os"
	"sort"
	"strings"
)

type NuclMsaJob struct {
	Size     int      `json:"size" validate:"required"`
	Database []string `json:"database"`
	Mode     string   `json:"mode" validate:"required"`
	query    string
}

func (r NuclMsaJob) Hash() Id {
	h := sha256.New224()
	h.Write([]byte(r.query))
	h.Write([]byte(r.Mode))

	sort.Strings(r.Database)

	for _, value := range r.Database {
		h.Write([]byte(value))
	}

	bs := h.Sum(nil)
	return Id(base64.URLEncoding.WithPadding(base64.NoPadding).EncodeToString(bs))
}

func (r NuclMsaJob) Rank() float64 {
	return float64(r.Size * max(len(r.Database), 1))
}

func (r NuclMsaJob) WriteFasta(path string) error {
	err := os.WriteFile(path, []byte(r.query), 0644)
	if err != nil {
		return err
	}
	return nil
}

func NewNuclMsaJobRequest(query string, dbs []string, validDbs []Params, mode string, resultPath string, email string) (JobRequest, error) {
	job := MsaJob{
		max(strings.Count(query, ">"), 1),
		dbs,
		mode,
		query,
	}

	request := JobRequest{
		job.Hash(),
		StatusPending,
		JobNuclMsa,
		job,
		email,
	}

	ids := make([]string, len(validDbs))
	for i, item := range validDbs {
		ids[i] = item.Path
	}

	for _, item := range job.Database {
		idx := isIn(item, ids)
		if idx == -1 {
			return request, errors.New("selected databases are not valid")
		}
	}

	return request, nil
}
