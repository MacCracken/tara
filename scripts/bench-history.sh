#!/bin/bash
set -euo pipefail

SAMPLES="${BENCH_SAMPLES:-100}"

echo "Running benchmarks (samples=$SAMPLES)..."
cargo bench --bench benchmarks -- --sample-size "$SAMPLES" --output-format bencher 2>/dev/null | tee bench-output.txt

# Append to CSV history
DATE=$(date -u +%Y-%m-%dT%H:%M:%SZ)
VERSION=$(cat VERSION | tr -d '[:space:]')
while IFS= read -r line; do
    echo "$DATE,$VERSION,$line" >> bench-history.csv
done < bench-output.txt

# Generate BENCHMARKS.md
{
    echo "# Benchmarks"
    echo ""
    echo "Version: $VERSION"
    echo "Date: $DATE"
    echo ""
    echo '```'
    cat bench-output.txt
    echo '```'
} > BENCHMARKS.md

rm -f bench-output.txt
echo "Done. Results in bench-history.csv and BENCHMARKS.md"
