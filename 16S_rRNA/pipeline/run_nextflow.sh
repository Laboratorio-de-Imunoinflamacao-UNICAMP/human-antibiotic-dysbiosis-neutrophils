#!/bin/bash

nextflow run nf-core/ampliseq \
	-params-file params_run.yaml \
	-profile docker \
	-r 2.12.0 \
	-resume \