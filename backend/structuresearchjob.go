package main

import (
	"crypto/sha256"
	"encoding/base64"
	"errors"
	"os"
	"sort"
	"strings"
)

type StructureSearchJob struct {
	Size            int      `json:"size" validate:"required"`
	Database        []string `json:"database" validate:"required"`
	Mode            string   `json:"mode" validate:"required,mode=3di tmalign 3diaa;print3di"`
	IterativeSearch bool     `json:"iterativesearch"`
	TaxFilter       string   `json:"taxfilter"`
	query           string
}

func (r StructureSearchJob) Hash() Id {
	h := sha256.New224()
	h.Write(([]byte)(JobStructureSearch))
	h.Write([]byte(r.query))
	h.Write([]byte(r.Mode))
	if r.IterativeSearch {
		h.Write([]byte("iterative"))
	}
	if r.TaxFilter != "" {
		h.Write([]byte(r.TaxFilter))
	}

	sort.Strings(r.Database)

	for _, value := range r.Database {
		h.Write([]byte(value))
	}

	bs := h.Sum(nil)
	return Id(base64.URLEncoding.WithPadding(base64.NoPadding).EncodeToString(bs))
}

func (r StructureSearchJob) Rank() float64 {
	modeFactor := 1
	if r.Mode == "tmalign" {
		modeFactor = 32
	}
	return float64(r.Size * max(len(r.Database), 1) * modeFactor)
}

func (r StructureSearchJob) WritePDB(path string) error {
	err := os.WriteFile(path, []byte(r.query), 0644)
	if err != nil {
		return err
	}
	return nil
}

func NewStructureSearchJobRequest(query string, dbs []string, validDbs []Params, mode string, resultPath string, email string, iterativeSearch bool, taxfilter string) (JobRequest, error) {
	job := StructureSearchJob{
		max(strings.Count(query, "HEADER"), 1),
		dbs,
		mode,
		iterativeSearch,
		taxfilter,
		query,
	}

	request := JobRequest{
		job.Hash(),
		StatusPending,
		JobStructureSearch,
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

	if !validTaxonFilter(taxfilter) {
		return request, errors.New("invalid taxon filter")
	}

	return request, nil
}
