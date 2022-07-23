.PHONY: build test

DOCARGS :=
JULIA ?= julia
mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

export HASTEST

build:
	cd $(realpath $(mkfile_path)..) && \
	$(JULIA) --startup-file=no --optimize=0 --project=docs/ \
		-e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()' && \
	$(JULIA) --startup-file=no --optimize=0 --project=docs/ \
		docs/make.jl ${DOCARGS}

test: DOCARGS+= --test
test: build
	@true

fix: DOCARGS+= --fix
fix: build
	@true
