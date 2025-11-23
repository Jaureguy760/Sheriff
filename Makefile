.PHONY: plots perf-suite

perf-suite:
	@echo "Running perf suite on superb slice..."
	python benchmarks/run_perf_suite.py \
		--dataset superb=benchmarks/large_superb.bam \
		--whitelist benchmarks/large_superb_whitelist.txt \
		--max-cells 30 \
		--max-genes 200 \
		--max-reads 100000 \
		--log benchmarks/run_perf_suite_superb.log

plots:
	@echo "Generating perf suite plots..."
	python benchmarks/plot_perf.py
